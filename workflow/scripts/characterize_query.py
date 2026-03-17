import csv
import shutil
import subprocess
from pathlib import Path


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


def parse_metrics_from_manifest(path: Path):
    metrics = {}
    with path.open() as handle:
        next(handle)
        for line in handle:
            field, value = line.rstrip("\n").split("\t", 1)
            metrics[field] = value
    return metrics


def count_gff_features(path: Path):
    if not path.exists():
        return 0
    count = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[2] == "CDS":
                count += 1
    return count


def run_prodigal(query_fasta: Path, proteins_out: Path, cds_out: Path, gff_out: Path):
    prodigal = shutil.which("prodigal")
    if not prodigal:
        raise RuntimeError("prodigal_not_found")
    cmd = [
        prodigal,
        "-i",
        str(query_fasta),
        "-a",
        str(proteins_out),
        "-d",
        str(cds_out),
        "-f",
        "gff",
        "-o",
        str(gff_out),
        "-p",
        "meta",
        "-q",
    ]
    subprocess.run(cmd, check=True)


def run_pyrodigal(query_fasta: Path, proteins_out: Path, cds_out: Path, gff_out: Path):
    try:
        import pyrodigal
    except ImportError as exc:
        raise RuntimeError(
            "No reusable query proteins were configured and neither 'prodigal' nor 'pyrodigal' is available."
        ) from exc

    with query_fasta.open() as handle:
        records = list(parse_fasta(query_fasta))

    finder = pyrodigal.GeneFinder(meta=True)
    with proteins_out.open("w") as proteins_handle, cds_out.open("w") as cds_handle, gff_out.open("w") as gff_handle:
        wrote_header = False
        for sequence_id, sequence in records:
            genes = finder.find_genes(sequence)
            genes.write_translations(
                proteins_handle,
                sequence_id=sequence_id.split()[0],
                include_stop=True,
                full_id=True,
            )
            genes.write_genes(
                cds_handle,
                sequence_id=sequence_id.split()[0],
                full_id=True,
            )
            genes.write_gff(
                gff_handle,
                sequence_id=sequence_id.split()[0],
                header=not wrote_header,
                full_id=True,
            )
            wrote_header = True


query_fasta = Path(snakemake.input.query)
manifest_path = Path(snakemake.input.manifest)
proteins_out = Path(snakemake.output.proteins)
cds_out = Path(snakemake.output.cds)
gff_out = Path(snakemake.output.gff)
metrics_out = Path(snakemake.output.metrics)

for path in [proteins_out.parent, cds_out.parent, gff_out.parent, metrics_out.parent]:
    path.mkdir(parents=True, exist_ok=True)

config_inputs = snakemake.config.get("inputs", {})
existing_proteins = Path(config_inputs["existing_query_proteins"]) if config_inputs.get("existing_query_proteins") else None
existing_cds = Path(config_inputs["existing_query_cds"]) if config_inputs.get("existing_query_cds") else None
existing_gff = Path(config_inputs["existing_query_gff"]) if config_inputs.get("existing_query_gff") else None
annotation_table = Path(config_inputs["existing_query_annotations"]) if config_inputs.get("existing_query_annotations") else None

source_mode = "reused_existing"
if existing_proteins and existing_proteins.exists():
    shutil.copyfile(existing_proteins, proteins_out)
    if existing_cds and existing_cds.exists():
        shutil.copyfile(existing_cds, cds_out)
    else:
        cds_out.write_text("")
    if existing_gff and existing_gff.exists():
        shutil.copyfile(existing_gff, gff_out)
    else:
        gff_out.write_text("")
else:
    try:
        source_mode = "predicted_with_prodigal"
        run_prodigal(query_fasta, proteins_out, cds_out, gff_out)
    except RuntimeError as exc:
        if str(exc) != "prodigal_not_found":
            raise
        source_mode = "predicted_with_pyrodigal"
        run_pyrodigal(query_fasta, proteins_out, cds_out, gff_out)

proteins = list(parse_fasta(proteins_out))
protein_count = len(proteins)
aa_total = sum(len(seq.rstrip("*")) for _, seq in proteins)
manifest_metrics = parse_metrics_from_manifest(manifest_path)
annotation_rows = 0
if annotation_table and annotation_table.exists():
    with annotation_table.open() as handle:
        annotation_rows = max(sum(1 for _ in handle) - 1, 0)

gff_cds_count = count_gff_features(gff_out)

rows = [
    ("metric", "value"),
    ("query_fasta", str(query_fasta.resolve())),
    ("query_records", manifest_metrics.get("query_records", "")),
    ("query_total_bp", manifest_metrics.get("query_total_bp", "")),
    ("orf_source", source_mode),
    ("protein_fasta", str(proteins_out.resolve())),
    ("protein_count", str(protein_count)),
    ("protein_total_aa", str(aa_total)),
    ("cds_fasta", str(cds_out.resolve())),
    ("cds_records", str(sum(1 for _ in parse_fasta(cds_out)) if cds_out.exists() else 0)),
    ("gff_file", str(gff_out.resolve())),
    ("gff_cds_features", str(gff_cds_count)),
    ("annotation_table", str(annotation_table.resolve()) if annotation_table and annotation_table.exists() else ""),
    ("annotation_rows", str(annotation_rows)),
]

with metrics_out.open("w") as handle:
    for key, value in rows:
        handle.write(f"{key}\t{value}\n")
