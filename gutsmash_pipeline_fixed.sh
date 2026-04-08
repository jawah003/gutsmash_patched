#!/bin/bash
# Fix 2: Removed `-i` (interactive shell flag) from shebang — not appropriate for batch jobs.
# Original: #!/bin/bash -i

# =============================================================================
# SLURM directives
# =============================================================================
#SBATCH --job-name=gutsmash_array
# Fix 7: Use full NFS paths for SLURM log files.
#SBATCH --output=/nfs/jjawahar/mi_gutsmash/logs/gutsmash_%A_%a.out
#SBATCH --error=/nfs/jjawahar/mi_gutsmash/logs/gutsmash_%A_%a.err
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
# Fix 6: Use the correct miniforge3 conda.sh path.
source ~/miniforge3/etc/profile.d/conda.sh
# Fix 5: Activate conda env by full path so the script works regardless of
#        whether the env name is registered in the current shell.
conda activate /nfs/jjawahar/miniforge3/envs/gutsmash_pipeline_v2/
set -u   # restore strict unbound-variable check

# =============================================================================
# Configuration
# =============================================================================
BASE_DIR="/nfs/jjawahar/mi_gutsmash"
GUTSMASH_DIR="${BASE_DIR}/gutsmash_patched"
MAGS_DIR="/nfs/abyrd/GNE/FDB1485/20_mags/nopurify"
RESULTS_DIR="${BASE_DIR}/results"

# Fix 3: Updated READS_LIST to the correct NFS path.
# Original: READS_LIST="${BASE_DIR}/forward_reads_list.txt"
# Can also copy over later
READS_LIST="/nfs/jjawahar/humann3_pipeline/2025-4-10_metabolomics_samples/fastq_file_list.txt"

# Fix 1: Replace generic antiSMASH rule IDs with the exact gutSMASH rule IDs
# sourced from antismash/detection/gut_hmm_detection/cluster_rules/strict.txt
# and relaxed.txt. These are used for exact set membership matching in Phase 2.
TARGET_PATHWAYS=(
    # --- strict.txt rules ---
    hydroxybenzoate2phenol pdu EUT_pathway TMA p-cresol
    Arginine2_Hcarbonate Arginine2putrescine acetate2butyrate
    Putrescine2spermidine proline2aminovalerate Leucine_reduction
    gallic_acid_met bai_operon AAA_reductive_branch porA
    PFOR_II_pathway Lysine_degradation glutamate2butyric
    caffeate_respiration carnitine_degradaion_caiTABCDE
    aminobutyrate2Butyrate succinate2propionate acrylate2propionate
    Threonine2propionate Pyruvate2acetate-formate Glycine_reductase
    Glycine_cleavage histidine2glutamate_hutHGIU_operon
    Oxidative_glycerol Acetyl-CoA_pathway Fumarate2succinate
    Indoleacetate2scatole Phenylacetate2toluene
    Hydroxy-L-proline2proline Sulfate2sulfide Rnf_complex
    Molybdopterin_dependent_oxidoreductase Nitrate_reductase
    Ech_complex Formate_dehydrogenase Respiratory_glycerol
    NADH_dehydrogenase_I Bilirubin_reductase
    Anaerobic_sulfite_redutase Phenylpyruvate_ferredoxin_oxidoreductase
    putative_2-oxoglutarate_ferredoxin_oxidoreductase
    Sulfate2PAPS PAPS2sulfide Adenylylsulfate_reductase
    r-butyrobetaine2TMA taurine2sulfite Alkanesulfonate2sulfite
    sulfoquinovose-EMP_pathway sulfoquinovose-ED_pathway
    xanthine_dehydrogenase uric_acid2SCFA
    # --- relaxed.txt rules ---
    Flavoenzyme_AA_peptides_catabolism Flavoenzyme_sugar_catabolism
    Flavoenzyme_lipids_catabolism OD_lactate_related
    OD_eut_pdu_related OD_AA_metabolism OD_fatty_acids
    OD_aldehydes_related OD_unknown TPP_fatty_acids TPP_AA_metabolism
    GR_AA_metabolism GR_eut-pdu-related GR_fatty_acids
    OD_GR_eut_related OD_GR_unassigned fatty_acids-unassigned
    Others_HGD_unassigned
)
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

# Fix 2: Create numbered subdirectories within scratch to match the reference
#        run_full_pipeline_array.sh structure for each pipeline phase.
WORK_DIR="${SCRATCH_DIR}/${SAMPLE_ID}"
DIR_GUTSMASH="${WORK_DIR}/1_gutsmash"
DIR_EXTRACT="${WORK_DIR}/2_cluster_extraction"
DIR_CDHIT="${WORK_DIR}/3_CD-HIT-EST"
DIR_SALMON="${WORK_DIR}/4_salmon_quant"
DIR_RPKM="${WORK_DIR}/5_abundance_RPKM"
mkdir -p "${DIR_GUTSMASH}" "${DIR_EXTRACT}" "${DIR_CDHIT}" "${DIR_SALMON}" "${DIR_RPKM}"

# =============================================================================
# Phase 1 — Run gutSMASH on all MAGs for this sample
# =============================================================================

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
    MAG_OUT_DIR="${DIR_GUTSMASH}/${MAG_BASENAME}"
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

EXTRACTED_FASTA="${DIR_EXTRACT}/${SAMPLE_ID}_clusters.fa"

# Fix 3: Pass the correct directory env vars to the Python heredoc.
export DIR_GUTSMASH EXTRACTED_FASTA

python3 - <<'PYTHON_SCRIPT'
# Fix 3 (Python extraction):
#   - Use DIR_GUTSMASH (numbered dir) instead of old SAMPLE_RESULTS_DIR.
#   - Feature type changed to "protocluster" — that is where gutSMASH stores
#     rule-matched clusters with product qualifiers (confirmed in run_full_pipeline.sh).
#   - Exact set membership matching: `p in targets` replaces substring matching.
#   - FASTA headers use "__" as delimiter (Salmon-safe; no "|" characters).
#   - Generates clusters_summary.csv alongside the FASTA.
import os
import sys
import json
import re
import glob
import csv

results_dir  = os.environ["DIR_GUTSMASH"]
output_fasta = os.environ["EXTRACTED_FASTA"]
targets      = set(os.environ.get("TARGET_PATHWAYS", "").split())

if not targets:
    raise RuntimeError("TARGET_PATHWAYS environment variable is empty or unset — aborting extraction")

extract_dir  = os.path.dirname(output_fasta)
summary_csv  = os.path.join(extract_dir, "clusters_summary.csv")

LOCATION_RE = re.compile(r"\[(\d+):(\d+)\]")

def parse_location(loc_str):
    """Return (start, end) integers from a location string like '[100:1200](+)'."""
    m = LOCATION_RE.search(loc_str)
    if not m:
        return None, None
    return int(m.group(1)), int(m.group(2))

written = 0
with open(output_fasta, "w") as out_fh, \
     open(summary_csv, "w", newline="") as csv_fh:

    writer = csv.writer(csv_fh)
    writer.writerow(["mag_id", "record_id", "product", "start", "end", "length"])

    for json_file in glob.glob(os.path.join(results_dir, "**", "*.json"), recursive=True):
        # mag_dir is the name of the gutSMASH output folder for this MAG
        mag_dir = os.path.basename(os.path.dirname(json_file))
        try:
            with open(json_file) as jf:
                data = json.load(jf)
        except (json.JSONDecodeError, OSError) as exc:
            print(f"WARNING: Could not parse {json_file}: {exc}", file=sys.stderr)
            continue

        for record in data.get("records", []):
            # Sequence is stored at the record level under record["seq"]["data"]
            record_seq = record.get("seq", {}).get("data", "")
            record_id  = record.get("id", "unknown")

            for feat in record.get("features", []):
                # Fix 3: Check for "protocluster" type — gutSMASH stores
                # rule-matched clusters here (not "region"/"subregion"/"cand_cluster").
                if feat.get("type", "") != "protocluster":
                    continue

                qualifiers = feat.get("qualifiers", {})
                product_list = qualifiers.get("product", [])
                products = product_list if isinstance(product_list, list) else [product_list]

                # Fix 3: Exact set membership — no substring/fuzzy matching.
                if not any(p in targets for p in products):
                    continue

                loc_str = feat.get("location", "")
                start, end = parse_location(loc_str)
                if start is None:
                    print(f"WARNING: Cannot parse location '{loc_str}' in {json_file}", file=sys.stderr)
                    continue

                cluster_seq = record_seq[start:end] if record_seq else "NNNN"
                if not cluster_seq:
                    cluster_seq = "NNNN"

                product = products[0] if products else "unknown"
                length  = end - start

                # Fix 4: Use "__" delimiter (Salmon-safe; avoids "|" parsing issues).
                header = f">{mag_dir}__{record_id}__{product}__{start}__{end}"
                out_fh.write(f"{header}\n{cluster_seq}\n")

                writer.writerow([mag_dir, record_id, product, start, end, length])
                written += 1

print(f"Extracted {written} cluster sequence(s) to {output_fasta}")
print(f"Summary written to {summary_csv}")
PYTHON_SCRIPT

# =============================================================================
# Phase 3 — Dereplication with CD-HIT-EST
# =============================================================================
echo "[$(date '+%F %T')] Phase 3: Dereplicating with CD-HIT-EST"

DEREP_FASTA="${DIR_CDHIT}/${SAMPLE_ID}_clusters_derep.fa"

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

SALMON_INDEX="${DIR_SALMON}/${SAMPLE_ID}_salmon_index"
SALMON_QUANT="${DIR_SALMON}/${SAMPLE_ID}_salmon_quant"

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

RPKM_CSV="${DIR_RPKM}/${SAMPLE_ID}_rpkm.csv"
INLINE_R_SCRIPT="${DIR_RPKM}/rpkm_normalize.R"

# Fix 5 (inlined R heredoc): Write the R normalisation script inline so the
# pipeline is fully self-contained (no external rpkm_normalize.R required).
# Args: 1=quant_file  2=out_file  3=total_library_reads
cat > "${INLINE_R_SCRIPT}" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript rpkm_normalize.R <quant.sf> <out.csv> <total_library_reads>")
}
quant_file          <- args[1]
out_file            <- args[2]
total_library_reads <- as.numeric(args[3])

quant <- read.table(quant_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# RPKM = (NumReads * 1e9) / (EffectiveLength * total_library_reads)
# Using total_library_reads (full library size) as denominator — NOT total_mapped —
# so values are comparable across samples regardless of pathway representation.
quant$RPKM <- (quant$NumReads * 1e9) / (quant$EffectiveLength * total_library_reads)

write.csv(quant[, c("Name", "NumReads", "TPM", "RPKM")], out_file, row.names = FALSE)
cat("RPKM normalisation complete:", out_file, "\n")
RSCRIPT

# Pass TOTAL_READS as the third argument so the R script uses it as the
# library-size denominator instead of total_mapped.
Rscript "${INLINE_R_SCRIPT}" \
    "${SALMON_QUANT}/quant.sf" \
    "${RPKM_CSV}" \
    "${TOTAL_READS}"

# =============================================================================
# Copy final results back to NFS (Fix 8: only finals leave local scratch)
# =============================================================================
echo "[$(date '+%F %T')] Copying results to NFS"

FINAL_DIR="${RESULTS_DIR}/${SAMPLE_ID}"
mkdir -p "${FINAL_DIR}"

cp "${SALMON_QUANT}/quant.sf"            "${FINAL_DIR}/${SAMPLE_ID}_quant.sf"
cp "${RPKM_CSV}"                         "${FINAL_DIR}/${SAMPLE_ID}_rpkm.csv"

# Copy clusters_summary.csv if it was generated (may be absent if no clusters matched)
SUMMARY_CSV="${DIR_EXTRACT}/clusters_summary.csv"
if [[ -f "${SUMMARY_CSV}" ]]; then
    cp "${SUMMARY_CSV}" "${FINAL_DIR}/clusters_summary.csv"
fi

echo "[$(date '+%F %T')] Done: ${SAMPLE_ID}"
