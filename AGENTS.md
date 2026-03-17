# AGENTS

## Purpose

This repository contains a reusable phage pangenome pipeline. The default user path is a single query FASTA plus automatic reference discovery.

## Default workflow

1. Put the query genome at `input/query/query.fasta`.
2. Run Snakemake with `workflow/Snakefile`.
3. Use `config/config.yaml` for the default remote-discovery mode.

## Supported discovery modes

- `blastn_remote`: default lightweight mode
- `local_blast_db`: optional faster/heavier mode using a prepared local BLAST database
- `configured`: use a prebuilt local cohort and GenBank file

## Editing rules

- Preserve the current config schema unless intentionally versioning a breaking change.
- Keep output paths under `results/<project_name>/`.
- Do not commit generated results, caches, virtual environments, or downloaded databases.
- Prefer extending `workflow/scripts/` and keeping rule files thin.
- Keep the default config portable; do not hardcode machine-specific project data into `config/config.yaml`.

## Versioning

- Use semantic tags such as `v0.1.0`, `v0.2.0`.
- Record user-visible changes in `CHANGELOG.md`.
- Treat config-schema changes and output-schema changes as version-relevant.
