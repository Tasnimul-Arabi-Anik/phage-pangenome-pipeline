# Generic Handover Rules for Any Phage Pangenome

## What to hand to another analyst

Always provide:

- query genome FASTA
- query protein FASTA
- candidate hit table from reference discovery
- final accession list for each cohort
- reference GenBank file used for the run
- table of excluded references with reasons
- orthogroup summary table
- presence/absence matrix
- query gene classification table
- report source

## What to tell them explicitly

- what lineage the pangenome is intended to represent
- whether the run is `close`, `expanded`, or `custom`
- whether near-duplicate references were removed
- whether zero-protein or incomplete GenBank records were excluded
- what orthology thresholds were used
- which genes are biologically important in the query

## Common phage-specific failure modes

- choosing references by host name instead of sequence similarity
- mixing distant phage groups in one pangenome
- trusting broken GenBank CDS annotations without QC
- treating singletons as biologically important before checking annotation quality
- ignoring partial ORFs at contig edges
- overinterpreting lysis proteins without topology or domain checks
