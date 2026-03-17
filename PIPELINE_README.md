# Snakemake Starter Pipeline

This is a starter scaffold for a generic phage comparative genomics and pangenome workflow.

## What it is

- A directory layout and rule structure for a reusable pipeline
- A config-driven `Snakefile`
- Rule stubs for each major stage of the workflow
- Example conda environment files

## What it is not

- It is not yet a fully implemented production workflow
- The early workflow stages are now script-backed, but later stages still need full scientific implementations

## Layout

- `workflow/Snakefile`: top-level workflow entry point
- `workflow/rules/`: one rule file per analysis stage
- `config/config.yaml`: main runtime configuration
- `envs/`: example conda environments
- `input/`: query and optional user-supplied references

## Default mode

The default config in [config.yaml](/home/anik/phage_project/pangenome/config/config.yaml) is now the validated `single FASTA + remote discovery` path.

Minimal user workflow:

1. Put the query FASTA at `input/query/query.fasta`
2. Run:

```bash
XDG_CACHE_HOME=$PWD/.cache ./.snakemake-venv/bin/snakemake -s workflow/Snakefile --cores 4
```

This default mode will:

- predict query ORFs with `prodigal` or `pyrodigal`
- run remote `blastn` against NCBI
- select candidate phage accessions by hit ranking and genome length filters
- download GenBank and FASTA records with NCBI `efetch`
- filter records lacking translated CDS features
- infer RBH orthogroups
- generate the standard heatmap and summary tables

## Alternative modes

The workflow currently supports three practical modes:

- `Configured local cohort`: set `inputs.existing_references_gb` and `inputs.existing_reference_accessions`
- `Single FASTA + remote discovery`: leave those reference inputs empty and set `reference_selection.discovery_method: blastn_remote`
- `Single FASTA + local BLAST DB`: leave those reference inputs empty, set `reference_selection.discovery_method: local_blast_db`, and provide `reference_selection.local_blast_db`

Use [local_cohort_example.yaml](/home/anik/phage_project/pangenome/config/local_cohort_example.yaml) as the starting point for the local-cohort mode.
Use [local_blast_db_example.yaml](/home/anik/phage_project/pangenome/config/local_blast_db_example.yaml) as the starting point for the optional local-database mode.

Example local-cohort launch:

```bash
XDG_CACHE_HOME=$PWD/.cache ./.snakemake-venv/bin/snakemake \
  -s workflow/Snakefile \
  --configfile config/local_cohort_example.yaml \
  --cores 4
```

Example local-BLAST-DB launch:

```bash
XDG_CACHE_HOME=$PWD/.cache ./.snakemake-venv/bin/snakemake \
  -s workflow/Snakefile \
  --configfile config/local_blast_db_example.yaml \
  --cores 4
```

## Validated status

Validated end-to-end:

- single FASTA plus remote NCBI discovery
- local configured cohort

Known caveat:

- `blastn -remote` can queue for several minutes before writing output
- `local_blast_db` avoids that queueing but requires a prepared local BLAST database

## Current implementation status

- Implemented:
  - `00_input.smk`: validates the configured query FASTA and writes an input manifest
  - `01_query.smk`: reuses precomputed query proteins/CDS/GFF or falls back to `prodigal` or `pyrodigal`
  - `02_discovery.smk`: supports configured cohorts, remote `blastn`, or local BLAST-database discovery
  - `03_retrieval.smk`: stages a local GenBank cohort or retrieves references from NCBI `efetch`
  - `04_qc.smk`: filters out reference records lacking translated CDS features and derives close/expanded cohorts
  - `05_proteins.smk`: extracts translated reference proteins and combines them with the query proteome
  - `06_orthology.smk`: infers RBH orthogroups and writes pangenome summary tables
  - `07_interpretation.smk`: classifies query genes into core, accessory, and singleton groups
  - `08_plots.smk`: generates the presence/absence heatmap in PNG and TIFF
  - `09_features.smk`: enriches query annotations from orthogroup consensus labels and writes a query-centered feature follow-up note
  - `10_report.smk`: builds manuscript-style Markdown and DOCX reports using the enriched annotation table

## Remaining work

1. Expand the feature-followup stage into modular `DeepTMHMM`, `SignalP`, and domain-search rules.
2. Enrich the report stage with more journal-style formatting and optional citations.
3. Add a more robust discovery backend if you want to avoid `blastn -remote` queueing.

## Important design rule

Keep biological decision logic in versioned code and config. Do not rely on memory or report prose alone to preserve reference-selection rules, exclusion rules, or orthology thresholds.
