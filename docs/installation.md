# Installation

## Overview

This pipeline is distributed as a Git repository plus a Snakemake workflow. It is not installed with `pip install phage-pangenome-pipeline`.

Typical use pattern:

1. clone the repository
2. create an environment
3. install Snakemake and required tools
4. place a query FASTA into `input/query/query.fasta`
5. run the workflow

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

## Alternative: Python virtual environment

```bash
git clone https://github.com/Tasnimul-Arabi-Anik/phage-pangenome-pipeline.git
cd phage-pangenome-pipeline

python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install snakemake pyrodigal pandas matplotlib seaborn
```

Then install the missing Python packages:

```bash
pip install pillow pyyaml
```

You still need external tools separately:

- `blastn`, `blastp`, `makeblastdb`
- `pandoc`

### Ubuntu or Debian example for external tools

```bash
sudo apt update
sudo apt install ncbi-blast+ pandoc
```

This route is acceptable, but conda or mamba is still preferred for the full workflow environment because it is easier to reproduce.

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

## First run

Place a FASTA file at `input/query/query.fasta` and run:

```bash
XDG_CACHE_HOME=$PWD/.cache snakemake -s workflow/Snakefile -n --cores 1
```

Then run the actual workflow:

```bash
XDG_CACHE_HOME=$PWD/.cache snakemake -s workflow/Snakefile --cores 4
```

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

## Troubleshooting

### `snakemake: command not found`

Your environment is not activated. Activate the `conda`, `mamba`, or `venv` environment first.

### `blastn: command not found`

BLAST+ is not installed or not available in `PATH`.

### `pandoc: command not found`

Pandoc is not installed or not available in `PATH`. The pipeline can still produce Markdown, but DOCX export requires Pandoc.

### Remote BLAST seems slow

This is normal. `blastn -remote` can wait several minutes before output appears.
