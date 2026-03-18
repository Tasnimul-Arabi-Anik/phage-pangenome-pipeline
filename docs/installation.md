# Installation

## Overview

This pipeline is distributed as a Git repository plus a Snakemake workflow. It is not installed with `pip install phage-pangenome-pipeline`.

Typical use pattern:

1. clone the repository
2. create an environment
3. install Snakemake and required tools
4. place a query FASTA into `input/query/query.fasta`
5. run a dry-run
6. run the workflow

Important distinction:

- Python packages such as `snakemake`, `pyrodigal`, and `pandas` can be installed with `pip`
- external command-line tools such as `blastn`, `blastp`, `makeblastdb`, and `pandoc` are not installed with `pip`
- for the `venv` route, install BLAST+ and Pandoc separately with your system package manager
- for the `conda` or `mamba` route, install them in the same environment

## Minimal requirements

- Python `3.12+`
- `snakemake`
- NCBI BLAST+ in `PATH`
- `pandoc` in `PATH` for DOCX report export
- Python packages required by the helper scripts:
  - `pyrodigal`
  - `pandas`
  - `matplotlib`
  - `seaborn`
  - `pillow`
  - `pyyaml`

## Recommended: mamba or conda

This is the best default for most bioinformatics users.

### Option A: one full environment

```bash
git clone https://github.com/Tasnimul-Arabi-Anik/phage-pangenome-pipeline.git
cd phage-pangenome-pipeline

mamba env create -f envs/full.yaml
mamba activate phage-pangenome-full
```

The `envs/full.yaml` environment includes:

- `snakemake`
- `pyrodigal`
- `pandas`
- `matplotlib`
- `seaborn`
- `pillow`
- `pyyaml`
- `blast`
- `pandoc`

If you use `conda` instead of `mamba`:

```bash
conda env create -f envs/full.yaml
conda activate phage-pangenome-full
```

### Option B: manual conda or mamba install

```bash
git clone https://github.com/Tasnimul-Arabi-Anik/phage-pangenome-pipeline.git
cd phage-pangenome-pipeline

mamba create -n phage-pangenome python=3.12 snakemake pyrodigal pandas matplotlib seaborn pillow pyyaml blast pandoc -c conda-forge -c bioconda
mamba activate phage-pangenome
```

This manual conda or mamba command is equivalent to the environment file and includes both Python dependencies and the required non-Python tools.

## Alternative: Python virtual environment

Use this route if you prefer standard Python tooling and are willing to install BLAST+ and Pandoc separately.

```bash
git clone https://github.com/Tasnimul-Arabi-Anik/phage-pangenome-pipeline.git
cd phage-pangenome-pipeline

python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install snakemake pyrodigal pandas matplotlib seaborn pillow pyyaml
```

What this installs:

- `snakemake`
- `pyrodigal`
- `pandas`
- `matplotlib`
- `seaborn`
- `pillow`
- `pyyaml`

What this does not install:

- `blastn`
- `blastp`
- `makeblastdb`
- `pandoc`

You still need those external tools separately.

### Ubuntu or Debian example for external tools

```bash
sudo apt update
sudo apt install ncbi-blast+ pandoc
```

If you want to use the optional Pfam domain-annotation stage, also install HMMER:

```bash
sudo apt install hmmer
```

This route is acceptable, but conda or mamba is still preferred for the full workflow environment because it is easier to reproduce.

### Practical fresh-clone `venv` example

This is the exact style of setup used in a fresh clone validation run:

```bash
git clone https://github.com/Tasnimul-Arabi-Anik/phage-pangenome-pipeline.git
cd phage-pangenome-pipeline

python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install snakemake pyrodigal pandas matplotlib seaborn pillow pyyaml

sudo apt update
sudo apt install ncbi-blast+ pandoc

cp /path/to/your/phage.fasta input/query/query.fasta
XDG_CACHE_HOME=$PWD/.cache snakemake -s workflow/Snakefile -n --cores 1
XDG_CACHE_HOME=$PWD/.cache snakemake -s workflow/Snakefile --cores 4
```

## Verify installation

Run:

```bash
python --version
snakemake --version
blastn -version
blastp -version
makeblastdb -version
pandoc --version | head -n 1
```

Expected result:

- all commands should be found
- Snakemake should print a version number
- BLAST+ tools should print version info
- Pandoc should print a version line

You can also confirm the Python helper modules explicitly:

```bash
python -c "import pyrodigal, pandas, matplotlib, seaborn, PIL, yaml; print('python dependencies OK')"
```

## First run

Place a FASTA file at `input/query/query.fasta` and run:

```bash
XDG_CACHE_HOME=$PWD/.cache snakemake -s workflow/Snakefile -n --cores 1
```

Then run the actual workflow:

```bash
XDG_CACHE_HOME=$PWD/.cache snakemake -s workflow/Snakefile --cores 4
```

If you used a local `venv`, make sure it is activated first:

```bash
source .venv/bin/activate
```

If you used `conda` or `mamba`, make sure that environment is activated first.

For the bundled smoke-test dataset:

```bash
XDG_CACHE_HOME=$PWD/.cache snakemake -s workflow/Snakefile -n --configfile tests/config/ci_smoke_test.yaml --cores 1
```

## Notes on BLAST backends

- `blastn_remote`
  Requires internet access and may queue at NCBI.
- `local_blast_db`
  Requires a prepared local BLAST database and more disk space, but avoids remote queueing.
- `configured`
  Requires a local accession list and a local multi-record GenBank file.

## Optional external BLASTP annotation

If you want stronger protein annotation than orthogroup consensus alone, provide a protein FASTA database and enable the optional external BLASTP stage in `config/config.yaml`:

```yaml
annotation:
  blastp: true
  protein_db_fasta: /path/to/annotation_database.faa
  protein_db_metadata: /path/to/annotation_metadata.tsv
  max_target_seqs: 5
  max_evalue: 1e-5
  min_identity: 30
  min_query_coverage: 50
```

This stage uses:

- `makeblastdb` to build a temporary protein BLAST database
- `blastp` to search the query proteins against that database
- optional metadata mapping to turn accession-only subject IDs into real product and module labels

Outputs:

- `features/query_blastp_hits.tsv`
- `features/query_blastp_summary.tsv`
- improved propagated labels in `features/query_gene_annotations.tsv`

Recommended metadata columns:

- `seq_id`
- `protein_id`
- `product`
- `module`

## Optional Pfam hmmscan annotation

If you want domain-level support for proteins that remain hypothetical after BLASTP and orthogroup consensus, enable the optional Pfam stage:

```yaml
annotation:
  pfam_hmmscan: true
  pfam_db: /path/to/Pfam-A.hmm
  pfam_max_evalue: 1e-3
```

This stage uses:

- `hmmscan` from HMMER
- a user-managed Pfam-format HMM database

Outputs:

- `features/query_hmmscan_hits.tsv`
- `features/query_hmmscan_summary.tsv`

This stage is optional and disabled by default because the Pfam database can require substantial disk space.

## Troubleshooting

### `snakemake: command not found`

Your environment is not activated. Activate the `conda`, `mamba`, or `venv` environment first.

### `blastn: command not found`

BLAST+ is not installed or not available in `PATH`.

If you used a `venv`, this is expected until you install BLAST+ separately with `apt`, `conda`, or another package manager.

### `pandoc: command not found`

Pandoc is not installed or not available in `PATH`. The pipeline can still produce Markdown, but DOCX export requires Pandoc.

If you used a `venv`, install it separately because `pip` does not provide the system `pandoc` executable used here.

### `hmmscan: command not found`

HMMER is not installed or not available in `PATH`. This matters only if you enable:

- `annotation.pfam_hmmscan: true`

For Ubuntu or Debian, install it with:

```bash
sudo apt install hmmer
```

### Query characterization fails because `prodigal` or `pyrodigal` is missing

The default single-FASTA mode requires one of these:

- `prodigal` available in `PATH`, or
- `pyrodigal` installed in the active Python environment

The recommended fix is to install `pyrodigal` in the same environment as Snakemake.

### Remote BLAST seems slow

This is normal. `blastn -remote` can wait several minutes before output appears.
