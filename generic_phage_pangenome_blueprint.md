# Generic Phage Pangenome Pipeline Blueprint

## Goal

Build a reusable pipeline for comparative genomics and pangenomics of a single query phage against a curated cohort of related reference phages. The pipeline should be usable across different phage systems, not only Klebsiella phages.

## Recommendation

Build this in three layers:

1. `Pipeline`
   The executable workflow that performs retrieval, filtering, clustering, plotting, and reporting.
2. `Method guide`
   The biological decision rules that tell users how to choose references and interpret results.
3. `Optional Codex skill`
   A thin wrapper for guided execution and report drafting after the pipeline is stable.

The pipeline is the core asset. A skill should sit on top of it, not replace it.

## Workflow manager

Preferred choices:

- `Snakemake`
  Best if you want a transparent, Python-friendly workflow with easy debugging.
- `Nextflow`
  Best if you want stronger portability, profiles, and containerized execution.

For a research lab workflow with moderate complexity, `Snakemake` is usually the better first implementation.

## Biological design principles

### 1. Define the pangenome scope before running anything

Each run should declare one of the following modes:

- `close`
  Immediate lineage pangenome using the closest dereplicated references.
- `expanded`
  Broader lineage pangenome using a larger BLAST-ranked, dereplicated cohort.
- `custom`
  User-supplied cohort.

Do not mix distant phage groups simply because they all infect the same host genus.

### 2. Prefer lineage-driven reference selection

Reference selection should start from the query itself:

- search the query against a nucleotide database
- rank hits by similarity and coverage
- restrict by genome quality and genome size
- remove obvious outliers
- dereplicate near-identical genomes

Do not build the cohort from keyword search alone.

### 3. Use a protein-based pangenome

For phages, protein-level clustering is usually more informative and more robust than nucleotide-level clustering because:

- gene order can be conserved while nucleotide divergence is substantial
- ORF boundaries can differ slightly
- lysis, tail, and hypothetical genes often evolve quickly

### 4. Track annotation provenance explicitly

For every query gene and every orthogroup, preserve:

- original ORF ID
- protein sequence source
- annotation source
- evidence type
- reference match
- orthogroup assignment

This becomes critical during manuscript writing.

## Proposed pipeline stages

### Stage 00. Input validation

Inputs:

- query genome FASTA
- optional user reference list
- optional user reference GenBank file

Checks:

- FASTA has one or more valid nucleotide records
- sequence lengths are plausible for phages
- filenames and sample names are normalized

Outputs:

- validated input manifest

### Stage 01. Query characterization

Tasks:

- compute genome length and GC content
- predict CDS with a phage-appropriate gene caller
- optionally predict tRNAs
- create query protein FASTA, CDS FASTA, and GFF

Tools:

- `Prodigal` or `Pyrodigal`
- `tRNAscan-SE` optionally

Outputs:

- `query_proteins.faa`
- `query_cds.fna`
- `query.gff`
- `query_metrics.tsv`

### Stage 02. Reference discovery

Tasks:

- search query genome against public nucleotide database
- rank hits by identity, coverage, and hit length
- keep only likely full genomes

Tools:

- remote `BLASTN`
- optional local cached hit tables
- optional preconfigured local accession lists for reproducible reruns

Outputs:

- `blast_hits.tsv`
- `candidate_accessions.tsv`

### Stage 03. Reference retrieval

Tasks:

- download GenBank and FASTA for candidate references
- or stage a prebuilt local multi-record GenBank cohort
- record accession, title, genome length, source

Outputs:

- `references.gb`
- `references.fna`
- `reference_metadata.tsv`

### Stage 04. Reference QC and dereplication

Tasks:

- remove genomes without usable CDS translations
- remove incomplete or suspicious records
- optionally dereplicate near-identical references
- split into `close` and `expanded` cohorts

Critical rules:

- exclude zero-protein records
- exclude records with broken annotations if protein pangenome is the goal
- document why references were removed

Outputs:

- `close_refs.txt`
- `expanded_refs.txt`
- `dropped_references.tsv`
- `cohort_metadata.tsv`

### Stage 05. Protein extraction and harmonization

Tasks:

- extract reference protein translations from GenBank
- harmonize IDs and metadata
- merge with query proteins

Outputs:

- `combined_proteins.faa`
- `protein_manifest.tsv`

### Stage 06. Orthology inference

Tasks:

- run all-vs-all protein comparison
- cluster orthologs into orthogroups
- classify orthogroups as core, accessory, singleton

Possible engines:

- `MMseqs2` preferred for scale and speed
- `BLASTP` acceptable for small to moderate cohorts
- optional `OrthoFinder` for a more formal orthology framework

Pragmatic note:

- a validated starter implementation can use reciprocal-best-hit `BLASTP` first
- switch to `MMseqs2` later if cohort sizes become large

Default thresholds should be configurable:

- minimum percent identity
- minimum alignment coverage
- maximum E-value

Outputs:

- `orthogroups.tsv`
- `presence_absence.tsv`
- `summary.tsv`

### Stage 07. Query-centered interpretation

Tasks:

- map each query protein to orthogroup category
- identify query singletons
- flag variable structural, lysis, and host-recognition genes

Outputs:

- `query_gene_orthogroup_classification.tsv`
- `query_singletons.tsv`
- `query_accessory_highlight.tsv`

### Stage 08. Plotting

Recommended standard figures:

- orthogroup presence/absence heatmap
- core/accessory/singleton summary plot
- optional lysis-cassette or feature-specific schematic

Outputs:

- `pangenome_presence_absence_heatmap.png`
- `pangenome_presence_absence_heatmap.tiff`
- other figure files as needed

### Stage 09. Feature-specific follow-up analyses

This stage is optional and should be modular.

Examples:

- DeepTMHMM for lysis proteins
- SignalP for secretory or SAR-like regions
- InterProScan or HHpred for difficult hypothetical proteins
- targeted BLASTP for singleton genes

Outputs:

- feature-specific result folders

### Stage 10. Report generation

Tasks:

- build journal-ready tables
- draft methods and results summaries
- export markdown and docx reports

Outputs:

- `report.md`
- `report.docx`
- `important_outputs/`

## Proposed directory structure

```text
phage_pangenome/
  workflow/
    Snakefile
    rules/
      00_input.smk
      01_query.smk
      02_discovery.smk
      03_retrieval.smk
      04_qc.smk
      05_proteins.smk
      06_orthology.smk
      07_interpretation.smk
      08_plots.smk
      09_features.smk
      10_report.smk
  config/
    config.yaml
    tools.yaml
  envs/
    base.yaml
    blast.yaml
    orthology.yaml
    plotting.yaml
  scripts/
    fetch_references.py
    extract_genbank_proteins.py
    cluster_orthogroups.py
    plot_heatmap.py
    build_report.py
  input/
    query/
    user_references/
  results/
    close/
    expanded/
  docs/
    methods_guide.md
    cohort_selection_guide.md
    interpretation_guide.md
```

## Configuration model

The pipeline should be driven by a small YAML config such as:

```yaml
project_name: example_phage
query_fasta: input/query/query.fasta
mode: expanded

reference_selection:
  max_candidates: 100
  expanded_target: 50
  close_target: 15
  min_genome_length: 15000
  max_genome_length: 500000
  require_translated_cds: true
  dereplicate: true

orthology:
  method: mmseqs2
  min_identity: 30
  min_query_coverage: 70
  min_subject_coverage: 70
  max_evalue: 1e-5

plots:
  make_heatmap: true
  image_dpi: 400

feature_modules:
  run_deeptmhmm: true
  run_signalp: false
```

## Minimum viable implementation

If you want to build this incrementally, implement in this order:

1. input validation
2. query ORF calling
3. BLAST-based reference discovery
4. GenBank retrieval
5. protein extraction
6. orthogroup clustering
7. presence/absence table
8. heatmap
9. query singleton/accessory summary
10. report export

Do not start with figure polish or a skill.

## What should be documented for every run

Each analysis should preserve:

- query accession or sample ID
- date of reference retrieval
- database queried
- candidate hit table
- excluded references and reasons
- dereplication logic
- orthology thresholds
- software versions
- final cohort membership

This is what makes the pangenome transferable to another analyst.

## When a Codex skill becomes useful

A skill is worth building after the pipeline exists and stabilizes.

That skill should:

- ask for the query FASTA
- ask for `close` or `expanded` mode
- run the pipeline
- summarize core/accessory/singleton results
- draft report text
- warn about common cohort-selection mistakes

It should not contain the scientific logic only in prose. The logic should remain in code and documented configuration.

## Final recommendation

For generic phage work:

- build a `Snakemake` or `Nextflow` pipeline first
- keep cohort selection and interpretation rules in versioned markdown docs
- add a skill only after the pipeline is reproducible

That gives you a system another person can actually rerun, audit, and extend.
