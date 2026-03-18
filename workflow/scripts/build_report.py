import argparse
import csv
from collections import Counter
from pathlib import Path


def read_kv_tsv(path: Path):
    data = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            data[row["metric"]] = row["value"]
    return data


def read_rows(path: Path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def md_table(headers, rows):
    if not rows:
        return ""
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(header, "")) for header in headers) + " |")
    return "\n".join(lines)


def parse_args():
    parser = argparse.ArgumentParser(description="Build a manuscript-style pangenome report.")
    parser.add_argument("--summary", required=True)
    parser.add_argument("--query-table", required=True)
    parser.add_argument("--annotation-summary", required=True)
    parser.add_argument("--feature-note", required=True)
    parser.add_argument("--heatmap", required=True)
    parser.add_argument("--output-md", required=True)
    parser.add_argument("--project-name", required=True)
    parser.add_argument("--mode", required=True)
    return parser.parse_args()


args = parse_args()

summary = read_kv_tsv(Path(args.summary))
query_rows = read_rows(Path(args.query_table))
annotation_summary = read_kv_tsv(Path(args.annotation_summary))
feature_note = Path(args.feature_note).read_text().strip()

report_path = Path(args.output_md)
report_path.parent.mkdir(parents=True, exist_ok=True)

project_name = args.project_name
mode = summary.get("mode", args.mode)
genomes = int(summary["genomes"])
reference_genomes = int(summary["reference_genomes"])
orthogroups_total = int(summary["orthogroups_total"])
core = int(summary["core_orthogroups"])
accessory = int(summary["accessory_orthogroups"])
singleton = int(summary["singleton_orthogroups"])
query_total = int(summary["query_proteins_total"])
query_core = int(summary["query_core_orthogroups"])
query_accessory = int(summary["query_accessory_orthogroups"])
query_singleton = int(summary["query_singleton_orthogroups"])

module_counter = Counter((row.get("preferred_module") or row.get("module") or "unassigned") for row in query_rows)
singletons = [row for row in query_rows if row["category"] == "singleton"]
accessory_rows = [row for row in query_rows if row["category"] == "accessory"]
accessory_rows = sorted(
    accessory_rows,
    key=lambda row: (-int(row["n_genomes"]), row.get("preferred_module", row.get("module", "")), row["gene_id"]),
)

summary_table = [
    {"Metric": "Total genomes", "Value": genomes},
    {"Metric": "Reference genomes", "Value": reference_genomes},
    {"Metric": "Total orthogroups", "Value": orthogroups_total},
    {"Metric": "Core orthogroups", "Value": core},
    {"Metric": "Accessory orthogroups", "Value": accessory},
    {"Metric": "Singleton orthogroups", "Value": singleton},
    {"Metric": "Query proteins", "Value": query_total},
    {"Metric": "Query core assignments", "Value": query_core},
    {"Metric": "Query accessory assignments", "Value": query_accessory},
    {"Metric": "Query singleton assignments", "Value": query_singleton},
]

top_accessory_table = []
for row in accessory_rows[:10]:
    top_accessory_table.append(
        {
            "Gene": row["gene_id"],
            "Module": row.get("preferred_module") or row["module"] or "",
            "Product": row.get("preferred_product") or row["product"] or row["consensus_product"] or "unannotated protein",
            "Genomes": row["n_genomes"],
        }
    )

singleton_table = []
for row in singletons:
    singleton_table.append(
        {
            "Gene": row["gene_id"],
            "Module": row.get("preferred_module") or row["module"] or "",
            "Product": row.get("preferred_product") or row["product"] or row["consensus_product"] or "unannotated protein",
        }
    )

unassigned_count = module_counter.get("unassigned", 0)
known_module_count = query_total - unassigned_count
core_fraction = (query_core / query_total) if query_total else 0.0
accessory_fraction = (query_accessory / query_total) if query_total else 0.0
singleton_fraction = (query_singleton / query_total) if query_total else 0.0

interpretation_lines = [
    "## Interpretation",
    "",
    f"The query phage is represented in a `{mode}` pangenome containing `{reference_genomes}` reference genomes and `{genomes}` genomes total. The query retains a largely conserved backbone, with `{query_core}` of `{query_total}` proteins ({core_fraction:.1%}) assigned to core orthogroups, while `{query_accessory}` proteins ({accessory_fraction:.1%}) fall into the accessory fraction and `{query_singleton}` protein ({singleton_fraction:.1%}) remains query-specific under the current clustering thresholds.",
    "",
]

if known_module_count == 0:
    interpretation_lines.extend(
        [
            "Because this run started from a raw FASTA without curated functional annotation, module labels and product names remain mostly generic. This is expected for a prediction-only run and does not indicate a pipeline error. Richer interpretation requires either curated annotation input or downstream domain-based annotation of the predicted proteins.",
            "",
        ]
    )
else:
    interpretation_lines.extend(
        [
            f"Functional module assignments were available for `{known_module_count}` query proteins, allowing a more specific biological interpretation of the conserved and variable gene fractions.",
            "",
        ]
    )

resolved_from_consensus = int(annotation_summary.get("resolved_from_consensus", "0"))
resolved_from_blastp = int(annotation_summary.get("annotation_source::external_blastp", "0"))
resolved_from_pfam = int(annotation_summary.get("annotation_source::pfam_hmmscan", "0"))
if resolved_from_consensus:
    interpretation_lines.extend(
        [
            f"Orthogroup consensus propagation improved annotation for `{resolved_from_consensus}` query proteins that would otherwise have remained generic or hypothetical in a FASTA-only run.",
            "",
        ]
    )
if resolved_from_blastp:
    interpretation_lines.extend(
        [
            f"External BLASTP annotation contributed preferred functional labels for `{resolved_from_blastp}` query proteins, providing a stronger annotation layer than pangenome consensus alone when a suitable protein database is configured.",
            "",
        ]
    )
if resolved_from_pfam:
    interpretation_lines.extend(
        [
            f"Pfam hmmscan contributed domain-based annotation for `{resolved_from_pfam}` query proteins, helping to classify proteins that remained weakly annotated after direct label propagation.",
            "",
        ]
    )

if singletons:
    singleton_products = ", ".join(
        row.get("preferred_product") or row["product"] or row["consensus_product"] or row["gene_id"]
        for row in singletons[:3]
    )
    interpretation_lines.extend(
        [
            f"The query-specific singleton set includes: {singleton_products}. These genes are the highest-priority targets for follow-up annotation, topology prediction, and homolog review because they are most likely to capture lineage-specific biology.",
            "",
        ]
    )

interpretation_lines.extend(
    [
        "## Conclusion",
        "",
        "This report should be read as a first-pass comparative pangenome summary. Core and accessory structure is already informative from a FASTA-only input, and orthogroup-consensus propagation provides a useful first annotation pass, but product-level biological interpretation improves substantially when the query is accompanied by curated annotation, domain predictions, or targeted follow-up analyses for singleton and accessory genes.",
    ]
)

methods_lines = [
    "## Methods",
    "",
    "The pipeline was run from a single query phage FASTA. Coding sequences were predicted with `prodigal` when available or with `pyrodigal` as the local fallback. Reference phages were identified either from a configured local cohort or by remote `blastn` discovery against NCBI nucleotide records, followed by GenBank retrieval with `efetch`. Records lacking translated CDS features were excluded before protein clustering.",
    "",
    f"Orthogroups were inferred with an all-vs-all reciprocal-best-hit `blastp` strategy using minimum identity `{summary['min_identity']}%`, minimum query coverage `{summary['min_query_coverage']}%`, minimum subject coverage `{summary['min_subject_coverage']}%`, and maximum E-value `{summary['max_evalue']}`. Orthogroups were classified as core, accessory, or singleton according to the number of genomes represented in each cluster.",
]
if resolved_from_blastp:
    methods_lines.extend(
        [
            "",
            "An optional external annotation step was enabled in which query proteins were searched against a user-supplied protein FASTA with `blastp`, and informative top hits were used to refine query product labels before report generation.",
        ]
    )
if resolved_from_pfam:
    methods_lines.extend(
        [
            "",
            "An optional Pfam annotation step was enabled in which query proteins were searched against a user-supplied HMM database with `hmmscan`, and informative domain hits were used as an additional annotation source.",
        ]
    )

results_lines = [
    "## Results",
    "",
    f"The `{project_name}` run was executed in `{mode}` mode and retained `{reference_genomes}` reference genomes plus the query genome (`{genomes}` genomes total). Across the final protein set, the pipeline identified `{orthogroups_total}` orthogroups, including `{core}` core groups, `{accessory}` accessory groups, and `{singleton}` singleton groups.",
    "",
    f"The query phage contributed `{query_total}` proteins to the pangenome. Of these, `{query_core}` mapped to core orthogroups, `{query_accessory}` mapped to accessory orthogroups, and `{query_singleton}` mapped to singleton orthogroups.",
    "",
    "### Summary statistics",
    "",
    md_table(["Metric", "Value"], summary_table),
    "",
    "### Query module distribution",
    "",
]

for module, count in sorted(module_counter.items(), key=lambda item: (-item[1], item[0])):
    results_lines.append(f"- `{module}`: {count}")

results_lines.extend(
    [
        "",
        "### Query accessory genes",
        "",
        md_table(["Gene", "Module", "Product", "Genomes"], top_accessory_table) or "No query accessory genes were detected.",
        "",
        "### Query singleton genes",
        "",
        md_table(["Gene", "Module", "Product"], singleton_table) or "No query singleton genes were detected.",
        "",
        "## Figures",
        "",
        "Figure 1. Orthogroup presence/absence heatmap generated from the final pangenome matrix.",
        "",
        f"![Pangenome heatmap]({Path(args.heatmap).resolve()})",
        "",
        "## Feature follow-up note",
        "",
        feature_note,
        "",
        *interpretation_lines,
    ]
)

report_text = "\n".join(
    [
        "# Comparative Pangenome Report",
        "",
        f"Project: `{project_name}`",
        "",
        *methods_lines,
        "",
        *results_lines,
        "",
    ]
)

report_path.write_text(report_text)
