# Phage Pangenome Pipeline

Reusable Snakemake pipeline for comparative genomics and pangenomics of a single query phage genome.

[![release](https://img.shields.io/github/v/tag/Tasnimul-Arabi-Anik/phage-pangenome-pipeline?label=release)](https://github.com/Tasnimul-Arabi-Anik/phage-pangenome-pipeline/releases)
[![license](https://img.shields.io/github/license/Tasnimul-Arabi-Anik/phage-pangenome-pipeline)](LICENSE)
[![workflow](https://img.shields.io/badge/workflow-Snakemake-039BE5)](https://snakemake.readthedocs.io/)
[![python](https://img.shields.io/badge/python-3.12%2B-3776AB)](https://www.python.org/)

## Validated modes

- single FASTA plus remote NCBI discovery
- configured local cohort
- optional local BLAST-database backend for heavier repeated use

## What it does

- starts from a single phage FASTA
- predicts CDS with `prodigal` or `pyrodigal`
- discovers related references from either:
  - remote NCBI `blastn`
  - a user-supplied local BLAST database
  - a configured local cohort
- retrieves and filters references
- builds a protein-based pangenome with reciprocal-best-hit `blastp`
- writes summary tables, a presence/absence heatmap, and a report in `md` and `docx`

## Quick start

1. Put your query genome at `input/query/query.fasta`
2. Review `config/config.yaml`
3. Run:

```bash
XDG_CACHE_HOME=$PWD/.cache ./.snakemake-venv/bin/snakemake -s workflow/Snakefile --cores 4
```

Default mode:

- `single FASTA + remote discovery`
- config file: `config/config.yaml`

Alternative example configs:

- `config/local_cohort_example.yaml`
- `config/local_blast_db_example.yaml`

## Installation notes

This repository expects `snakemake`, `blast+`, and Python dependencies required by the helper scripts.

The workflow was developed and validated with:

- Python `3.12`
- Snakemake `9.x`
- BLAST+ available in `PATH`
- `pandoc` available in `PATH` for DOCX report export

## Discovery modes

- `blastn_remote`
  Default lightweight mode. Good for occasional runs, but remote BLAST can queue for several minutes.
- `local_blast_db`
  Faster and more reproducible for repeated use, but requires a prepared local BLAST database.
- `configured`
  Best when you already have a curated reference cohort and GenBank file.

## Expected input

Place the query genome at:

```text
input/query/query.fasta
```

## Output

Outputs are written under:

```text
results/<project_name>/
```

Key outputs:

- `orthology/summary.tsv`
- `orthology/presence_absence.tsv`
- `interpretation/query_gene_orthogroup_classification.tsv`
- `plots/pangenome_presence_absence_heatmap.png`
- `plots/pangenome_presence_absence_heatmap.tiff`
- `report/report.md`
- `report/report.docx`

## Repository guide

- Remote `blastn` can queue for several minutes.
- Local BLAST mode is heavier but avoids remote queueing.
- Large downloaded references and generated results should not be committed.
- `PIPELINE_README.md` contains the fuller workflow guide.
- `AGENTS.md` contains Codex-oriented repository instructions.

## Citation

If you use this repository, cite the tagged release and see [CITATION.cff](CITATION.cff).
