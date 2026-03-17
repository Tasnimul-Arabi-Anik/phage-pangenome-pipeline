import csv
from pathlib import Path


def parse_accessions(path: Path):
    accessions = []
    if not path.exists():
        return accessions
    with path.open() as handle:
        for raw_line in handle:
            value = raw_line.strip()
            if not value or value.startswith("#"):
                continue
            accessions.append(value.split("\t", 1)[0].split(".", 1)[0])
    return accessions


def parse_genbank_records(path: Path):
    text = path.read_text()
    for chunk in text.split("\n//"):
        chunk = chunk.strip()
        if chunk:
            yield chunk.splitlines()


def parse_genbank_qualifiers(lines, start_index):
    qualifiers = {}
    i = start_index
    current_key = None
    while i < len(lines):
        line = lines[i]
        if line.startswith("ORIGIN") or line.startswith("     gene") or line.startswith("     CDS") or (
            line.startswith("     ") and not line.startswith("                     ")
        ):
            break
        if line.startswith("                     /"):
            content = line.strip()[1:]
            if "=" in content:
                key, value = content.split("=", 1)
                value = value.strip()
                if value.startswith('"'):
                    value = value[1:]
                    parts = []
                    if value.endswith('"'):
                        qualifiers[key] = value[:-1]
                    else:
                        parts.append(value)
                        i += 1
                        while i < len(lines):
                            cont = lines[i].strip()
                            if cont.endswith('"'):
                                parts.append(cont[:-1])
                                break
                            parts.append(cont)
                            i += 1
                        qualifiers[key] = "".join(parts)
                else:
                    qualifiers[key] = value
                current_key = key
            else:
                qualifiers[content] = ""
                current_key = content
        elif current_key:
            qualifiers[current_key] += line.strip()
        i += 1
    return qualifiers, i


def collect_genome_stats(path: Path):
    genome_rows = []
    for lines in parse_genbank_records(path):
        version = ""
        definition = ""
        organism = ""
        feature_start = None
        for idx, line in enumerate(lines):
            if line.startswith("VERSION"):
                version = line.split()[1]
            elif line.startswith("DEFINITION"):
                definition = line[len("DEFINITION") :].strip().rstrip(".")
            elif line.startswith("  ORGANISM"):
                organism = line[len("  ORGANISM") :].strip()
            elif line.startswith("FEATURES"):
                feature_start = idx + 1
                break
        genome_id = version or definition or organism
        accession = genome_id.split(".", 1)[0]
        protein_count = 0
        i = feature_start or 0
        while i < len(lines):
            line = lines[i]
            if line.startswith("ORIGIN"):
                break
            if line.startswith("     CDS"):
                qualifiers, i = parse_genbank_qualifiers(lines, i + 1)
                translation = qualifiers.get("translation", "").replace(" ", "")
                if translation:
                    protein_count += 1
                continue
            i += 1
        genome_rows.append(
            {
                "accession": accession,
                "genome_id": genome_id,
                "genome_name": organism or definition or genome_id,
                "protein_count": protein_count,
            }
        )
    return genome_rows


reference_gb = Path(snakemake.input.gb)
metadata_path = Path(snakemake.input.metadata)
close_out = Path(snakemake.output.close_refs)
expanded_out = Path(snakemake.output.expanded_refs)
dropped_out = Path(snakemake.output.dropped)
cohort_meta_out = Path(snakemake.output.cohort_meta)
for path in [close_out.parent, expanded_out.parent, dropped_out.parent, cohort_meta_out.parent]:
    path.mkdir(parents=True, exist_ok=True)

requested_accessions = []
if metadata_path.exists():
    with metadata_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        requested_accessions = [row["accession"].split(".", 1)[0] for row in reader]

requested_set = set(requested_accessions)
genome_rows = collect_genome_stats(reference_gb)
usable_rows = [row for row in genome_rows if row["protein_count"] > 0 and (not requested_set or row["accession"] in requested_set)]
usable_rows.sort(key=lambda row: requested_accessions.index(row["accession"]) if row["accession"] in requested_set else len(requested_accessions))

close_target = int(snakemake.config["reference_selection"]["close_target"])
expanded_target = int(snakemake.config["reference_selection"]["expanded_target"])

close_rows = usable_rows[:close_target]
expanded_rows = usable_rows[:expanded_target]
mode = snakemake.config["mode"]
selected_mode_accessions = {row["accession"] for row in (close_rows if mode == "close" else expanded_rows)}

with close_out.open("w") as handle:
    for row in close_rows:
        handle.write(f"{row['accession']}\n")

with expanded_out.open("w") as handle:
    for row in expanded_rows:
        handle.write(f"{row['accession']}\n")

usable_set = {row["accession"] for row in usable_rows}
with dropped_out.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["accession", "reason"])
    for row in genome_rows:
        if requested_set and row["accession"] not in requested_set:
            writer.writerow([row["accession"], "not_in_requested_reference_list"])
        elif row["accession"] not in usable_set:
            writer.writerow([row["accession"], "no_translated_cds"])

with cohort_meta_out.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(
        [
            "accession",
            "genome_id",
            "genome_name",
            "protein_count",
            "selected_close",
            "selected_expanded",
            "selected_for_analysis",
        ]
    )
    for row in usable_rows:
        writer.writerow(
            [
                row["accession"],
                row["genome_id"],
                row["genome_name"],
                row["protein_count"],
                str(row["accession"] in {item["accession"] for item in close_rows}).lower(),
                str(row["accession"] in {item["accession"] for item in expanded_rows}).lower(),
                str(row["accession"] in selected_mode_accessions).lower(),
            ]
        )
