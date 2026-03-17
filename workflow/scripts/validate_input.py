from pathlib import Path


def parse_fasta_stats(path: Path):
    record_count = 0
    total_length = 0
    current = []
    had_header = False

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if had_header:
                    total_length += len("".join(current))
                    current = []
                had_header = True
                record_count += 1
                continue
            current.append(line)

    if had_header:
        total_length += len("".join(current))

    return record_count, total_length


query_path = Path(snakemake.input.query)
if not query_path.exists():
    raise FileNotFoundError(f"Configured query FASTA does not exist: {query_path}")

record_count, total_length = parse_fasta_stats(query_path)
if record_count == 0:
    raise ValueError(f"Query FASTA contains no records: {query_path}")

output_path = Path(snakemake.output.manifest)
output_path.parent.mkdir(parents=True, exist_ok=True)

mode = snakemake.config["mode"]
config_inputs = snakemake.config.get("inputs", {})
existing_refs = config_inputs.get("existing_references_gb", "")
reference_source = Path(existing_refs) if existing_refs else None
reference_status = "configured" if reference_source and reference_source.exists() else "missing"

rows = [
    ("field", "value"),
    ("project_name", snakemake.config["project_name"]),
    ("mode", mode),
    ("query_fasta", str(query_path.resolve())),
    ("query_records", str(record_count)),
    ("query_total_bp", str(total_length)),
    ("existing_references_gb", str(reference_source.resolve()) if reference_source and reference_source.exists() else str(reference_source or "")),
    ("existing_references_status", reference_status),
]

with output_path.open("w") as handle:
    for key, value in rows:
        handle.write(f"{key}\t{value}\n")
