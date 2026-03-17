import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Protein:
    seq_id: str
    genome_id: str
    genome_name: str
    protein_id: str
    locus_tag: str
    product: str
    sequence: str
    source: str
    module: str = ""


def parse_fasta(path: Path):
    header = None
    seq = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line)
    if header is not None:
        yield header, "".join(seq)


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


def load_query_annotations(path: Path):
    annotations = {}
    if not path.exists():
        return annotations
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            annotations[row["gene_id"]] = {
                "product": row.get("product", "").strip() or "hypothetical protein",
                "module": row.get("module", "").strip(),
            }
    return annotations


def load_query_proteins(path: Path, annotation_path: Path | None):
    proteins = []
    annotations = load_query_annotations(annotation_path) if annotation_path else {}
    genome_id = "query"
    genome_name = "query"
    for header, sequence in parse_fasta(path):
        gene_id = header.split()[0]
        annotation = annotations.get(gene_id, {})
        proteins.append(
            Protein(
                seq_id=f"{genome_id}|{gene_id}",
                genome_id=genome_id,
                genome_name=genome_name,
                protein_id=gene_id,
                locus_tag=gene_id,
                product=annotation.get("product", "hypothetical protein"),
                module=annotation.get("module", ""),
                sequence=sequence.rstrip("*"),
                source="query",
            )
        )
    return proteins


def load_accessions(path: Path):
    accessions = set()
    if not path.exists():
        return accessions
    with path.open() as handle:
        for raw_line in handle:
            value = raw_line.strip()
            if not value or value.startswith("#"):
                continue
            accessions.add(value.split("\t", 1)[0].split(".", 1)[0])
    return accessions


def load_reference_proteins(path: Path, selected_accessions: set[str] | None):
    proteins = []
    genome_stats = {}
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
        if not genome_id:
            continue
        accession_root = genome_id.split(".", 1)[0]
        if selected_accessions and accession_root not in selected_accessions:
            continue
        genome_name = organism or definition or genome_id
        genome_stats.setdefault(genome_id, {"genome_name": genome_name, "protein_count": 0})
        i = feature_start or 0
        while i < len(lines):
            line = lines[i]
            if line.startswith("ORIGIN"):
                break
            if line.startswith("     CDS"):
                qualifiers, i = parse_genbank_qualifiers(lines, i + 1)
                translation = qualifiers.get("translation", "").replace(" ", "")
                if translation:
                    protein_id = qualifiers.get("protein_id", qualifiers.get("locus_tag", "unknown"))
                    proteins.append(
                        Protein(
                            seq_id=f"{genome_id}|{protein_id}",
                            genome_id=genome_id,
                            genome_name=genome_name,
                            protein_id=protein_id,
                            locus_tag=qualifiers.get("locus_tag", ""),
                            product=qualifiers.get("product", "hypothetical protein"),
                            sequence=translation.rstrip("*"),
                            source="reference",
                        )
                    )
                    genome_stats[genome_id]["protein_count"] += 1
                continue
            i += 1
    return proteins, genome_stats


query_proteins = Path(snakemake.input.query_proteins)
reference_gb = Path(snakemake.input.gb)
cohort_meta_path = Path(snakemake.input.cohort_meta)
combined_out = Path(snakemake.output.combined)
manifest_out = Path(snakemake.output.manifest)

for path in [combined_out.parent, manifest_out.parent]:
    path.mkdir(parents=True, exist_ok=True)

config_inputs = snakemake.config.get("inputs", {})
annotation_path = Path(config_inputs["existing_query_annotations"]) if config_inputs.get("existing_query_annotations") else None

selected_accessions = None
if cohort_meta_path.exists():
    with cohort_meta_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        accessions = [row["accession"].split(".", 1)[0] for row in reader if row.get("selected_for_analysis", "").lower() == "true"]
    if accessions:
        selected_accessions = set(accessions)

query_records = load_query_proteins(query_proteins, annotation_path)
reference_records = []
genome_stats = {}
if reference_gb.exists() and reference_gb.stat().st_size > 0:
    reference_records, genome_stats = load_reference_proteins(reference_gb, selected_accessions)

all_records = query_records + reference_records
if not all_records:
    raise ValueError("No proteins were loaded for the combined protein FASTA.")

with combined_out.open("w") as fasta_handle:
    for record in all_records:
        fasta_handle.write(f">{record.seq_id}\n")
        for i in range(0, len(record.sequence), 80):
            fasta_handle.write(record.sequence[i : i + 80] + "\n")

with manifest_out.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(
        [
            "seq_id",
            "genome_id",
            "genome_name",
            "protein_id",
            "locus_tag",
            "product",
            "module",
            "source",
            "aa_length",
        ]
    )
    for record in all_records:
        writer.writerow(
            [
                record.seq_id,
                record.genome_id,
                record.genome_name,
                record.protein_id,
                record.locus_tag,
                record.product,
                record.module,
                record.source,
                len(record.sequence),
            ]
        )
