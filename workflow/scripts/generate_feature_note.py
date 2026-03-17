import csv
from collections import Counter
from pathlib import Path


query_table_path = Path(snakemake.input.query_table)
note_path = Path(snakemake.output.note)
note_path.parent.mkdir(parents=True, exist_ok=True)

rows = []
with query_table_path.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    rows = list(reader)

category_counts = Counter(row["category"] for row in rows)
module_counts = Counter(row["module"] or "unassigned" for row in rows)
singletons = [row for row in rows if row["category"] == "singleton"]
accessory = [row for row in rows if row["category"] == "accessory"]
accessory_sorted = sorted(
    accessory,
    key=lambda row: (-int(row["n_genomes"]), row["module"], row["gene_id"]),
)

lines = [
    "# Feature Follow-up",
    "",
    "## Query orthogroup profile",
    "",
    f"- Core query genes: `{category_counts.get('core', 0)}`",
    f"- Accessory query genes: `{category_counts.get('accessory', 0)}`",
    f"- Singleton query genes: `{category_counts.get('singleton', 0)}`",
    "",
    "## Module distribution",
    "",
]

for module, count in sorted(module_counts.items(), key=lambda item: (-item[1], item[0])):
    lines.append(f"- `{module}`: {count}")

lines.extend(["", "## Suggested follow-up targets", ""])

if singletons:
    for row in singletons:
        product = row["product"] or row["consensus_product"] or "unannotated protein"
        lines.append(
            f"- Singleton: `{row['gene_id']}` ({product}); prioritize topology/domain searches and homolog review."
        )
else:
    lines.append("- No query singleton genes were detected.")

if accessory_sorted:
    for row in accessory_sorted[:8]:
        product = row["product"] or row["consensus_product"] or "unannotated protein"
        lines.append(
            f"- Accessory: `{row['gene_id']}` ({product}); present in `{row['n_genomes']}` genomes."
        )

note_path.write_text("\n".join(lines) + "\n")
