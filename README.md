# Phage Pangenome Pipeline

Reusable Snakemake pipeline for comparative genomics and pangenomics of a single query phage genome.

[![release](https://img.shields.io/github/v/tag/Tasnimul-Arabi-Anik/phage-pangenome-pipeline?label=release)](https://github.com/Tasnimul-Arabi-Anik/phage-pangenome-pipeline/releases)
[![license](https://img.shields.io/github/license/Tasnimul-Arabi-Anik/phage-pangenome-pipeline)](LICENSE)
[![workflow](https://img.shields.io/badge/workflow-Snakemake-039BE5)](https://snakemake.readthedocs.io/)
[![python](https://img.shields.io/badge/python-3.12%2B-3776AB)](https://www.python.org/)
[![ci-smoke](https://github.com/Tasnimul-Arabi-Anik/phage-pangenome-pipeline/actions/workflows/ci-smoke.yml/badge.svg)](https://github.com/Tasnimul-Arabi-Anik/phage-pangenome-pipeline/actions/workflows/ci-smoke.yml)

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
- enriches query annotations from orthogroup consensus labels for FASTA-only runs
- optionally runs external `BLASTP` annotation against a user-supplied protein FASTA database
- optionally runs `hmmscan` against a user-supplied Pfam HMM database
- writes summary tables, a presence/absence heatmap, and a report in `md` and `docx`

## Workflow overview

See [docs/workflow_overview.md](docs/workflow_overview.md) for a compact architecture diagram and stage summary.

## Quick start

1. Clone the repository
2. Create an environment and install dependencies
3. Put your query genome at `input/query/query.fasta`
4. Review `config/config.yaml`
5. Run a dry-run
6. Run the full workflow

```bash
git clone https://github.com/Tasnimul-Arabi-Anik/phage-pangenome-pipeline.git
cd phage-pangenome-pipeline
```

Recommended environment setup with `mamba`:

```bash
mamba env create -f envs/full.yaml
mamba activate phage-pangenome-full
```

Alternative setup with a local Python virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install snakemake pyrodigal pandas matplotlib seaborn pillow pyyaml
```

For the `venv` route, install the non-Python tools separately. `blastn` and `pandoc`
do not come from `pip`:

```bash
sudo apt update
sudo apt install ncbi-blast+ pandoc
```

If you want to use the optional Pfam domain-annotation stage, also install HMMER:

```bash
sudo apt install hmmer
```

Dry-run first:

```bash
XDG_CACHE_HOME=$PWD/.cache snakemake -s workflow/Snakefile -n --cores 1
```

Then run the workflow:

```bash
XDG_CACHE_HOME=$PWD/.cache snakemake -s workflow/Snakefile --cores 4
```

Default mode:

- `single FASTA + remote discovery`
- config file: `config/config.yaml`

Alternative example configs:

- `config/local_cohort_example.yaml`
- `config/local_blast_db_example.yaml`

## Adjust reference cohort size

The default single-FASTA configuration is conservative and uses:

- `max_candidates: 50`
- `close_target: 10`
- `expanded_target: 20`

This means the pipeline may discover up to 50 candidate phages, but the default expanded analysis keeps only the top 20 references after filtering and QC.

To use more references, edit `config/config.yaml`:

```yaml
reference_selection:
  max_candidates: 80
  close_target: 15
  expanded_target: 40
```

Practical guidance:

- use `10-20` references for a tight close-relative analysis
- use `30-50` references for a broader lineage-level pangenome
- increase `max_candidates` when you increase `expanded_target`, so the pipeline has enough filtered candidates to choose from
- if `require_translated_cds: true`, some downloaded references may still be excluded before final selection

After editing the config, rerun:

```bash
XDG_CACHE_HOME=$PWD/.cache snakemake -s workflow/Snakefile --cores 4
```

## Installation notes

This repository expects `snakemake`, `blast+`, and Python dependencies required by the helper scripts.

The workflow was developed and validated with:

- Python `3.12`
- Snakemake `9.x`
- BLAST+ available in `PATH`
- `pandoc` available in `PATH` for DOCX report export

Recommended full environment file:

- `envs/full.yaml`

See [docs/installation.md](docs/installation.md) for the full installation guide, including:

- exact fresh-clone setup
- `mamba` or `conda`
- Python `venv` plus `pip`
- BLAST+ installation
- Pandoc installation
- verification commands
- first dry-run and first full run

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
- `features/query_blastp_hits.tsv`
- `features/query_blastp_summary.tsv`
- `features/query_hmmscan_hits.tsv`
- `features/query_hmmscan_summary.tsv`
- `features/query_gene_annotations.tsv`
- `features/annotation_summary.tsv`
- `plots/pangenome_presence_absence_heatmap.png`
- `plots/pangenome_presence_absence_heatmap.tiff`
- `report/report.md`
- `report/report.docx`

## Test dataset

A small bundled test genome is provided at:

```text
tests/data/klebsiella.fasta
```

Lightweight CI uses:

```text
tests/config/ci_smoke_test.yaml
```

## Repository guide

- Remote `blastn` can queue for several minutes.
- Local BLAST mode is heavier but avoids remote queueing.
- FASTA-only runs now receive a first-pass annotation boost from orthogroup consensus propagation.
- For stronger annotation, enable `annotation.blastp: true` and set `annotation.protein_db_fasta` to a protein FASTA database.
- If your BLAST database FASTA headers are accession-only, also provide `annotation.protein_db_metadata` with `seq_id` or `protein_id` to product mappings.
- For domain-level support, enable `annotation.pfam_hmmscan: true` and set `annotation.pfam_db` to a Pfam-format HMM database.
- External BLASTP annotation improves labels, but it still does not replace domain-based annotation such as InterProScan or HMMER.
- Large downloaded references and generated results should not be committed.
- `PIPELINE_README.md` contains the fuller workflow guide.
- `AGENTS.md` contains Codex-oriented repository instructions.

## Citation

If you use this repository, cite the tagged release and see [CITATION.cff](CITATION.cff).

Suggested citation:

```text
Anik TA. Phage Pangenome Pipeline. GitHub repository. Version v0.1.0. 2026.
```

BibTeX:

```bibtex
@software{anik_phage_pangenome_pipeline_2026,
  author  = {Anik, Tasnimul Arabi},
  title   = {Phage Pangenome Pipeline},
  year    = {2026},
  version = {v0.1.0},
  url     = {https://github.com/Tasnimul-Arabi-Anik/phage-pangenome-pipeline}
}
```
