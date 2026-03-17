# Installation

## Minimal requirements

- Python `3.12+`
- `snakemake`
- NCBI BLAST+ in `PATH`
- `pandoc` in `PATH` for DOCX report export

## Option 1: local virtual environment

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install snakemake pyrodigal pandas matplotlib seaborn
```

Install external tools separately:

- `blastn`, `blastp`, `makeblastdb`
- `pandoc`

## Option 2: conda or mamba

```bash
mamba create -n phage-pangenome python=3.12 snakemake pandas matplotlib seaborn pyrodigal -c conda-forge -c bioconda
mamba activate phage-pangenome
mamba install blast pandoc -c bioconda -c conda-forge
```

## Quick validation

Place a FASTA file at `input/query/query.fasta` and run:

```bash
XDG_CACHE_HOME=$PWD/.cache snakemake -s workflow/Snakefile -n --cores 1
```

For the bundled smoke-test dataset:

```bash
XDG_CACHE_HOME=$PWD/.cache snakemake -s workflow/Snakefile -n --configfile tests/config/ci_smoke_test.yaml --cores 1
```
