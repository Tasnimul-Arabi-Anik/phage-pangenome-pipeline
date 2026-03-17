# Changelog

## Unreleased

- Added an optional external `BLASTP` annotation stage using a user-supplied protein FASTA database.
- Added optional `annotation.protein_db_metadata` support so accession-only BLAST databases can still resolve to product and module labels.
- Added `features/query_blastp_hits.tsv` and `features/query_blastp_summary.tsv`.
- Updated report generation to distinguish annotation coming from query labels, orthogroup consensus, and external BLASTP.

## v0.2.0

- Added annotation enrichment for FASTA-only runs by propagating orthogroup consensus labels into query-gene outputs.
- Added `features/query_gene_annotations.tsv` and `features/annotation_summary.tsv`.
- Improved manuscript-style reports with explicit interpretation and conclusion sections.
- Improved feature follow-up summaries with annotation-source counts and enriched preferred labels.
- Documented how to increase close and expanded reference cohort sizes in the README.

## v0.1.0

- Initial public export of the phage pangenome pipeline.
- Added validated single-FASTA remote-discovery mode.
- Added optional local-cohort mode.
- Added optional local BLAST-database discovery mode.
- Added RBH `blastp` orthogroup inference, heatmap generation, and report generation.
