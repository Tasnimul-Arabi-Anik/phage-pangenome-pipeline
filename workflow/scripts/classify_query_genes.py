import csv
from pathlib import Path


presence_absence_path = Path(snakemake.input.presence_absence)
manifest_path = Path(snakemake.input.manifest)
query_table_path = Path(snakemake.output.query_table)
query_singletons_path = Path(snakemake.output.query_singletons)
query_table_path.parent.mkdir(parents=True, exist_ok=True)

query_manifest = {}
with manifest_path.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        if row["source"] == "query":
            query_manifest[row["protein_id"]] = row

query_rows = []
singleton_rows = []
with presence_absence_path.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        query_members = [value.strip() for value in row.get("query", "").split(",") if value.strip()]
        for gene_id in query_members:
            meta = query_manifest.get(gene_id, {})
            query_rows.append(
                {
                    "orthogroup": row["orthogroup"],
                    "category": row["category"],
                    "n_genomes": row["n_genomes"],
                    "gene_id": gene_id,
                    "module": meta.get("module", ""),
                    "product": meta.get("product", ""),
                    "consensus_product": row["consensus_product"],
                }
            )
            if row["category"] == "singleton":
                singleton_rows.append(
                    {
                        "gene_id": gene_id,
                        "orthogroup": row["orthogroup"],
                        "product": meta.get("product", row["consensus_product"]),
                    }
                )

query_rows.sort(key=lambda item: ({"core": 0, "accessory": 1, "singleton": 2}[item["category"]], -int(item["n_genomes"]), item["gene_id"]))

with query_table_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["orthogroup", "category", "n_genomes", "gene_id", "module", "product", "consensus_product"],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(query_rows)

with query_singletons_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["gene_id", "orthogroup", "product"], delimiter="\t")
    writer.writeheader()
    writer.writerows(singleton_rows)
