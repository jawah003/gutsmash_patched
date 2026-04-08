#!/bin/bash
# Fix 2: Removed `-i` (interactive shell flag) from shebang — not appropriate for batch jobs.
# Original: #!/bin/bash -i

# =============================================================================
# SLURM directives
# =============================================================================
#SBATCH --job-name=gutsmash_array
#SBATCH --output=logs/gutsmash_%A_%a.out
#SBATCH --error=logs/gutsmash_%A_%a.err
#SBATCH --time=24:00:00
# Fix 1: Changed --mem=200G → --mem=16G.
# Original --mem=200G requested 200 GB *per task*, limiting concurrency to 1.
# gutSMASH + Salmon peak ~8-10 GB sequentially; 16 GB gives comfortable headroom.
# At 16 GB per task, 12 simultaneous jobs = 192 GB (fits within the 200 GB pool).
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
# Fix 1 (cont.): Changed --array=1-100%10 → 1-2%12 for initial testing.
# NOTE: Change to --array=1-441%12 for full production run (441 samples, 12 concurrent).
#SBATCH --array=1-2%12

# Fix 2 (cont.): Add strict mode immediately after SBATCH directives.
# set -e  : exit on any command failure
# set -u  : treat unset variables as errors
# set -o pipefail : catch failures inside pipelines
# NOTE: set +u / set -u brackets around `conda activate` to avoid spurious
#       "unbound variable" errors from conda's internal scripts.
set -euo pipefail

# =============================================================================
# Environment setup
# =============================================================================
set +u   # temporarily relax unbound-variable check for conda activation
# shellcheck source=/dev/null
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gutsmash_env
set -u   # restore strict unbound-variable check

# =============================================================================
# Configuration
# =============================================================================
BASE_DIR="/nfs/jjawahar/gutsmash_pipeline"
GUTSMASH_DIR="${BASE_DIR}/gutsmash_patched"
MAGS_DIR="${BASE_DIR}/mags"
RESULTS_DIR="${BASE_DIR}/results"

# Fix 3: Updated READS_LIST to the correct NFS path.
# Original: READS_LIST="${BASE_DIR}/forward_reads_list.txt"
READS_LIST="/nfs/jjawahar/humann3_pipeline/2025-4-10_metabolomics_samples/fastq_file_list.txt"

TARGET_PATHWAYS=("RiPP" "NRPS" "PKS" "terpene" "saccharide")
# Export as a space-separated string so the Python heredoc can read it via
# os.environ["TARGET_PATHWAYS"] — bash arrays cannot be exported directly.
export TARGET_PATHWAYS="${TARGET_PATHWAYS[*]}"

# Fix 8: Use node-local scratch ($TMPDIR, or /tmp as fallback) for intermediate I/O.
# Running 12 concurrent jobs writing gutSMASH/CD-HIT/Salmon intermediates to NFS
# can saturate network filesystem bandwidth. Compute on local disk, copy finals back.
SCRATCH_DIR="${TMPDIR:-/tmp}/gutsmash_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${SCRATCH_DIR}"

# Ensure scratch is cleaned up on exit (normal, error, or signal).
# shellcheck disable=SC2064
trap "rm -rf '${SCRATCH_DIR}'" EXIT

# =============================================================================
# Resolve sample for this array task
# =============================================================================
FORWARD_READ=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${READS_LIST}")
# Fix 3: Derive reverse read from the forward read path using the documented pattern.
# Forward reads end in _clean.1.fastq.gz; reverse reads end in _clean.2.fastq.gz.
REVERSE_READ="${FORWARD_READ/.1.fastq.gz/.2.fastq.gz}"
SAMPLE_ID=$(basename "${FORWARD_READ}" "_clean.1.fastq.gz")

echo "[$(date '+%F %T')] Starting sample: ${SAMPLE_ID}"
echo "  Forward read : ${FORWARD_READ}"
echo "  Reverse read : ${REVERSE_READ}"
echo "  Scratch dir  : ${SCRATCH_DIR}"

# =============================================================================
# Phase 1 — Run gutSMASH on all MAGs for this sample
# =============================================================================
SAMPLE_RESULTS_DIR="${SCRATCH_DIR}/gutsmash_results"
mkdir -p "${SAMPLE_RESULTS_DIR}"

# Fix 7: Use nullglob so that a pattern with no matching files expands to nothing
#        instead of being kept as a literal string in the array.
# Original: SAMPLE_MAGS=(${MAGS_DIR}/${SAMPLE_ID}*.{fa,fna,fasta} 2>/dev/null || true)
# (The 2>/dev/null on glob expansion has no effect; || true applies to the assignment,
#  not the glob; unexpanded literals would enter the array.)
shopt -s nullglob
SAMPLE_MAGS=("${MAGS_DIR}/${SAMPLE_ID}"*.fa \
             "${MAGS_DIR}/${SAMPLE_ID}"*.fna \
             "${MAGS_DIR}/${SAMPLE_ID}"*.fasta)
shopt -u nullglob

if [[ ${#SAMPLE_MAGS[@]} -eq 0 ]]; then
    echo "WARNING: No MAG files found for sample ${SAMPLE_ID}. Skipping." >&2
    exit 0
fi

echo "[$(date '+%F %T')] Phase 1: Running gutSMASH on ${#SAMPLE_MAGS[@]} MAG(s)"

for MAG in "${SAMPLE_MAGS[@]}"; do
    MAG_BASENAME=$(basename "${MAG%.*}")
    MAG_OUT_DIR="${SAMPLE_RESULTS_DIR}/${MAG_BASENAME}"
    mkdir -p "${MAG_OUT_DIR}"

    echo "  Processing MAG: ${MAG_BASENAME}"
    # Fix 4: Added --cpus "$SLURM_CPUS_PER_TASK" to gutSMASH invocation.
    # Without this flag, gutSMASH defaults to multiprocessing.cpu_count() (ALL node
    # cores), oversubscribing the node and violating SLURM's resource contract.
    # Source: antismash/config/args.py lines 466-470; flows to hmmsearch via
    # antismash/common/subprocessing/hmmsearch.py line 30.
    python3 "${GUTSMASH_DIR}/run_gutsmash.py" \
        --cpus "${SLURM_CPUS_PER_TASK}" \
        --output-dir "${MAG_OUT_DIR}" \
        --genefinding-tool prodigal \
        "${MAG}"
done

# =============================================================================
# Phase 2 — Extract target pathway cluster sequences from gutSMASH JSON output
# =============================================================================
echo "[$(date '+%F %T')] Phase 2: Extracting target cluster sequences"

EXTRACTED_FASTA="${SCRATCH_DIR}/${SAMPLE_ID}_clusters.fa"

python3 - <<'PYTHON_SCRIPT'
# Fix 5A: The original loop used `feature.get(...)` inside `for feat in features:`.
#         `feature` was never defined — this caused an immediate NameError.
# Fix 5B: The original used `feature.get("sequence", "NNNN")` to retrieve the
#         nucleotide sequence. In gutSMASH's JSON (antiSMASH serialisation),
#         sequences are stored at the *record* level under record["seq"]["data"],
#         NOT inside feature dicts. The fallback "NNNN" always fired, producing
#         an output FASTA full of NNNN strings.
#
# Correct approach:
#   1. Read the parent record's sequence from record["seq"]["data"].
#   2. Parse the feature's "location" string (format: "[start:end](strand)") to
#      obtain start/end coordinates.
#   3. Slice the record sequence to get the actual cluster nucleotides.
import os
import sys
import json
import re
import glob

results_dir  = os.environ["SAMPLE_RESULTS_DIR"]
output_fasta = os.environ["EXTRACTED_FASTA"]
targets      = set(os.environ.get("TARGET_PATHWAYS", "RiPP NRPS PKS terpene saccharide").split())

LOCATION_RE = re.compile(r"\[(\d+):(\d+)\]")

def parse_location(loc_str):
    """Return (start, end) integers from a location string like '[100:1200](+)'."""
    m = LOCATION_RE.search(loc_str)
    if not m:
        return None, None
    return int(m.group(1)), int(m.group(2))

written = 0
with open(output_fasta, "w") as out_fh:
    for json_file in glob.glob(os.path.join(results_dir, "**", "*.json"), recursive=True):
        try:
            with open(json_file) as jf:
                data = json.load(jf)
        except (json.JSONDecodeError, OSError) as exc:
            print(f"WARNING: Could not parse {json_file}: {exc}", file=sys.stderr)
            continue

        for record in data.get("records", []):
            # Fix 5B: sequence is at the record level, under record["seq"]["data"]
            record_seq = record.get("seq", {}).get("data", "")
            record_id  = record.get("id", "unknown")

            for feat in record.get("features", []):        # Fix 5A: was `feature`
                feat_type = feat.get("type", "")           # Fix 5A: was `feature.get`
                if feat_type not in ("region", "subregion", "cand_cluster"):
                    continue

                qualifiers = feat.get("qualifiers", {})
                # gutSMASH/antiSMASH stores the BGC type in the "product" qualifier
                product_list = qualifiers.get("product", [])
                products = product_list if isinstance(product_list, list) else [product_list]

                # Check whether any product matches a target pathway
                if not any(any(t.lower() in p.lower() for t in targets) for p in products):
                    continue

                loc_str = feat.get("location", "")
                start, end = parse_location(loc_str)
                if start is None:
                    print(f"WARNING: Cannot parse location '{loc_str}' in {json_file}", file=sys.stderr)
                    continue

                # Fix 5B: slice the record-level sequence for the actual nucleotides
                cluster_seq = record_seq[start:end] if record_seq else "NNNN"
                if not cluster_seq:
                    cluster_seq = "NNNN"

                product_str = "_".join(products) if products else "unknown"
                header = f">{record_id}|{feat_type}|{start}-{end}|{product_str}"
                out_fh.write(f"{header}\n{cluster_seq}\n")
                written += 1

print(f"Extracted {written} cluster sequence(s) to {output_fasta}")
PYTHON_SCRIPT

# =============================================================================
# Phase 3 — Dereplication with CD-HIT-EST
# =============================================================================
echo "[$(date '+%F %T')] Phase 3: Dereplicating with CD-HIT-EST"

DEREP_FASTA="${SCRATCH_DIR}/${SAMPLE_ID}_clusters_derep.fa"

cd-hit-est \
    -i "${EXTRACTED_FASTA}" \
    -o "${DEREP_FASTA}" \
    -c 0.95 \
    -T "${SLURM_CPUS_PER_TASK}" \
    -M 14000 \
    -d 0

# =============================================================================
# Phase 4 — Quantification with Salmon
# =============================================================================
echo "[$(date '+%F %T')] Phase 4: Quantifying with Salmon"

SALMON_INDEX="${SCRATCH_DIR}/${SAMPLE_ID}_salmon_index"
SALMON_QUANT="${SCRATCH_DIR}/${SAMPLE_ID}_salmon_quant"

salmon index \
    -t "${DEREP_FASTA}" \
    -i "${SALMON_INDEX}" \
    --threads "${SLURM_CPUS_PER_TASK}"

salmon quant \
    -i "${SALMON_INDEX}" \
    -l A \
    -1 "${FORWARD_READ}" \
    -2 "${REVERSE_READ}" \
    -p "${SLURM_CPUS_PER_TASK}" \
    --validateMappings \
    -o "${SALMON_QUANT}"

# =============================================================================
# Phase 5 — RPKM normalisation with R
# =============================================================================
echo "[$(date '+%F %T')] Phase 5: RPKM normalisation"

# Fix 6: Count total reads in the original FASTQ to use as the RPKM denominator.
# The original R script divided by total_mapped (reads mapped only to extracted
# pathway sequences), which is scientifically invalid: samples with low pathway
# representation get artificially inflated values.
# We must normalise by total *library* size (all reads in the input FASTQ).
#
# NOTE: `zcat | awk 'END{print NR/4}'` is slow for large files (must decompress
# the whole file) but is the most reliable method — `zgrep -c '^@'` can over-count
# because '@' also appears in FASTQ quality strings.
# If read counts are pre-computed (e.g., from a QC log or a companion .count file),
# replace this line with a simple `cat` to avoid the decompression cost.
echo "  Counting total reads in ${FORWARD_READ} ..."
TOTAL_READS=$(zcat "${FORWARD_READ}" | awk 'END{print NR/4}')
echo "  Total reads: ${TOTAL_READS}"

RPKM_CSV="${SCRATCH_DIR}/${SAMPLE_ID}_rpkm.csv"

# Pass TOTAL_READS as the third argument so the R script uses it as the
# library-size denominator (args[3]) instead of total_mapped.
Rscript "${BASE_DIR}/rpkm_normalize.R" \
    "${SALMON_QUANT}/quant.sf" \
    "${RPKM_CSV}" \
    "${TOTAL_READS}"

# =============================================================================
# Copy final results back to NFS (Fix 8: only finals leave local scratch)
# =============================================================================
echo "[$(date '+%F %T')] Copying results to NFS"

FINAL_DIR="${RESULTS_DIR}/${SAMPLE_ID}"
mkdir -p "${FINAL_DIR}"

cp "${SALMON_QUANT}/quant.sf" "${FINAL_DIR}/${SAMPLE_ID}_quant.sf"
cp "${RPKM_CSV}"              "${FINAL_DIR}/${SAMPLE_ID}_rpkm.csv"

echo "[$(date '+%F %T')] Done: ${SAMPLE_ID}"
