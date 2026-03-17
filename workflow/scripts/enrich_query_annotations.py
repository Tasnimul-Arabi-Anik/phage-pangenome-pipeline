import csv
from collections import Counter
from pathlib import Path


HYPOTHETICAL_TERMS = {
    "",
    "hypothetical protein",
    "conserved hypothetical protein",
    "putative protein",
    "uncharacterized protein",
    "unknown function",
}


MODULE_PATTERNS = [
    ("lysis", ["holin", "endolysin", "lysin", "spanin", "rz", "rzi", "amidase", "muramidase"]),
    ("packaging", ["terminase", "portal", "packaging"]),
    ("head", ["capsid", "head", "scaffold", "procapsid", "decoration"]),
    ("tail", ["tail", "fiber", "baseplate", "neck", "tape measure", "tube"]),
    ("dna_metabolism", ["polymerase", "primase", "helicase", "exonuclease", "recombinase", "methyltransferase", "kinase", "dna", "ssdna"]),
    ("transcription_regulation", ["transcriptional regulator", "regulator", "repressor", "antiterminator"]),
    ("membrane", ["membrane", "transmembrane", "signal peptide"]),
]


def is_hypothetical(value: str) -> bool:
    normalized = value.strip().lower()
    return normalized in HYPOTHETICAL_TERMS


def infer_module(product: str) -> str:
    value = product.strip().lower()
    for module, patterns in MODULE_PATTERNS:
        if any(pattern in value for pattern in patterns):
            return module
    return "hypothetical" if is_hypothetical(product) else "other"


query_table_path = Path(snakemake.input.query_table)
annotation_table_path = Path(snakemake.output.annotation_table)
summary_path = Path(snakemake.output.summary)
annotation_table_path.parent.mkdir(parents=True, exist_ok=True)

with query_table_path.open() as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

enriched_rows = []
annotation_source_counts = Counter()
module_counts = Counter()
resolved_from_consensus = 0

for row in rows:
    original_product = (row.get("product") or "").strip()
    consensus_product = (row.get("consensus_product") or "").strip()
    original_module = (row.get("module") or "").strip()

    if original_product and not is_hypothetical(original_product):
        preferred_product = original_product
        annotation_source = "query_annotation"
    elif consensus_product and not is_hypothetical(consensus_product):
        preferred_product = consensus_product
        annotation_source = "orthogroup_consensus"
        resolved_from_consensus += 1
    elif original_product:
        preferred_product = original_product
        annotation_source = "query_annotation"
    elif consensus_product:
        preferred_product = consensus_product
        annotation_source = "orthogroup_consensus"
    else:
        preferred_product = "hypothetical protein"
        annotation_source = "none"

    if original_module and original_module != "hypothetical":
        preferred_module = original_module
    else:
        preferred_module = infer_module(preferred_product)

    annotation_source_counts[annotation_source] += 1
    module_counts[preferred_module] += 1

    enriched = dict(row)
    enriched["preferred_product"] = preferred_product
    enriched["preferred_module"] = preferred_module
    enriched["annotation_source"] = annotation_source
    enriched_rows.append(enriched)

with annotation_table_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "orthogroup",
            "category",
            "n_genomes",
            "gene_id",
            "module",
            "product",
            "consensus_product",
            "preferred_module",
            "preferred_product",
            "annotation_source",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(enriched_rows)

with summary_path.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["metric", "value"])
    writer.writerow(["query_genes", len(enriched_rows)])
    writer.writerow(["resolved_from_consensus", resolved_from_consensus])
    for source, count in sorted(annotation_source_counts.items()):
        writer.writerow([f"annotation_source::{source}", count])
    for module, count in sorted(module_counts.items()):
        writer.writerow([f"preferred_module::{module}", count])
