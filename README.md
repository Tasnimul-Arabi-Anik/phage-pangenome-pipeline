# Phage Pangenome Pipeline

Reusable Snakemake pipeline for comparative genomics and pangenomics of a single query phage genome.

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

## Main entry point

```bash
XDG_CACHE_HOME=$PWD/.cache ./.snakemake-venv/bin/snakemake -s workflow/Snakefile --cores 4
```

Default config:

- `config/config.yaml`

Alternative example configs:

- `config/local_cohort_example.yaml`
- `config/local_blast_db_example.yaml`

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

## Notes

- Remote `blastn` can queue for several minutes.
- Local BLAST mode is heavier but avoids remote queueing.
- Large downloaded references and generated results should not be committed.

See `PIPELINE_README.md` for more operational detail.
