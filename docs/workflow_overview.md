# Workflow Overview

## High-level flow

```mermaid
flowchart TD
    A[Query phage FASTA] --> B[ORF prediction<br/>prodigal or pyrodigal]
    B --> C{Reference discovery mode}
    C --> D[Remote BLASTN]
    C --> E[Local BLAST database]
    C --> F[Configured local cohort]
    D --> G[Reference retrieval]
    E --> G
    F --> H[Stage local GenBank cohort]
    G --> I[Reference QC and filtering]
    H --> I
    I --> J[Protein extraction and harmonization]
    B --> J
    J --> K[All-vs-all BLASTP]
    K --> L[RBH orthogroup inference]
    L --> M[Query gene classification]
    B --> S[Optional external BLASTP annotation]
    B --> T[Optional Pfam hmmscan annotation]
    M --> R[Query annotation enrichment]
    S --> R
    T --> R
    L --> N[Presence/absence heatmap]
    L --> O[Summary tables]
    R --> P[Feature follow-up note]
    N --> Q[Markdown and DOCX report]
    O --> Q
    P --> Q
```

## Main outputs

- `orthology/summary.tsv`
- `orthology/presence_absence.tsv`
- `interpretation/query_gene_orthogroup_classification.tsv`
- `features/query_blastp_hits.tsv`
- `features/query_hmmscan_hits.tsv`
- `features/query_gene_annotations.tsv`
- `features/annotation_summary.tsv`
- `plots/pangenome_presence_absence_heatmap.png`
- `report/report.md`
- `report/report.docx`

## Supported discovery backends

- `blastn_remote`
- `local_blast_db`
- `configured`
