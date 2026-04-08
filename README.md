# gutSMASH (Modern HPC & Conda Fork)

> **Fork Notice:** This is a modernized fork of the original [gutSMASH repository](https://github.com/victoriapascal/gutsmash). It has been patched to support Python 3.11+, Biopython 1.85, and modern Conda/HPC environments without requiring root `apt-get` privileges. 

Anaerobic bacteria in the gut are responsible for the synthesis and transformation of diverse molecules involved in host-microbe and microbe-microbe interactions. The pathways for the production of these molecules are often physically clustered in the genome as metabolic gene clusters (MGCs).

gutSMASH is a tool developed to systematically evaluate the metabolic potential of these bacteria by predicting both known and novel anaerobic MGCs from the gut microbiome.

---

## 🛠️ Modernization & Compatibility Patches

The original gutSMASH release pinned several legacy dependencies (e.g., Python 3.6, Biopython 1.76, scikit-learn 0.19.0). To allow this pipeline to run on modern High-Performance Computing (HPC) clusters, the following patches were applied to the source code:

* **Biopython 1.85 Compatibility:** Updated legacy API calls. Removed deprecated `Bio.Alphabet` stubs, rewrote `Seq` API usages, and updated `UnknownSeq` handling across serialiser and record processing modules.
* **Muscle v5 Support:** Patched the `muscle` subprocessing module (`antismash/common/subprocessing/muscle.py`) with a case-insensitive version check to support `muscle=5.3` provided via Bioconda.
* **Scikit-Learn Update:** Confirmed compatibility with `scikit-learn=1.7.1` when running under `--minimal` mode. *(Note: If utilizing ML-based cluster scoring without `--minimal`, you may need to downgrade to `0.22.1` due to pickle compatibility).*
* **Diamond Subprocessing:** Patched `diamond.py` to skip `--threads/--tmpdir` flags during the `version` subcommand to comply with newer Diamond releases.
* **Dependency Trimming:** Omitted `MOODS-python` from the base environment as it is only required when passing `--cb-knownclusters` or `--enable-genefunctions`. 

## 📦 Installation (Conda / Micromamba)

Instead of relying on system-level `apt-get` installations, this fork uses an `environment.yml` to pull all necessary binaries (Prodigal, HMMER, Diamond, Muscle, etc.) via Bioconda.

1. **Clone the repository:**
   ```bash
   git clone [https://github.com/](https://github.com/)<YOUR_USERNAME>/gutsmash.git
   cd gutsmash
   ```

2. **Build and activate the environment:**
	```bash
	# Using conda
	conda env create -f environment.yml
	conda activate gutsmash_pipeline

	# Or using micromamba (recommended for speed)
	micromamba env create -f environment.yml
	micromamba activate gutsmash_pipeline
	```
3. **Download databases:**
	```text
	gutSMASH requires specific databases. Download them from the official gutSMASH page and place them in the correct directory path as detailed there.
	```



