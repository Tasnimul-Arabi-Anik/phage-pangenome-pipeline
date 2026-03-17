import csv
import shutil
import subprocess
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
            accessions.append(value.split("\t", 1)[0])
    return accessions


def run(cmd):
    subprocess.run(cmd, check=True)


def download_efetch(accessions, rettype: str, output_path: Path):
    chunk_size = int(snakemake.config["reference_selection"]["remote_fetch_chunk_size"])
    with output_path.open("w") as out_handle:
        for index in range(0, len(accessions), chunk_size):
            chunk = accessions[index : index + chunk_size]
            ids = ",".join(chunk)
            url = (
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                f"?db=nucleotide&id={ids}&rettype={rettype}&retmode=text"
            )
            result = subprocess.run(["curl", "-LfsS", url], check=True, capture_output=True, text=True)
            out_handle.write(result.stdout)


config_inputs = snakemake.config.get("inputs", {})
source_gb = config_inputs.get("existing_references_gb", "")
source_accessions_value = config_inputs.get("existing_reference_accessions", "")
source_accessions = Path(source_accessions_value) if source_accessions_value else Path("")
candidates = parse_accessions(Path(snakemake.input.candidates))

gb_out = Path(snakemake.output.gb)
fna_out = Path(snakemake.output.fna)
metadata_out = Path(snakemake.output.metadata)
for path in [gb_out.parent, fna_out.parent, metadata_out.parent]:
    path.mkdir(parents=True, exist_ok=True)

metadata_rows = []
if source_gb:
    source_gb_path = Path(source_gb)
    if not source_gb_path.exists():
        raise FileNotFoundError(f"Configured reference GenBank file does not exist: {source_gb_path}")
    shutil.copyfile(source_gb_path, gb_out)
    fna_out.write_text("")
    accessions = parse_accessions(source_accessions)
    for accession in accessions:
        metadata_rows.append([accession, "configured", str(source_gb_path.resolve()), str(source_accessions.resolve())])
else:
    if not candidates:
        raise ValueError("No reference accessions were available for remote retrieval.")
    download_efetch(candidates, "gbwithparts", gb_out)
    download_efetch(candidates, "fasta", fna_out)
    for accession in candidates:
        metadata_rows.append([accession, "downloaded", "NCBI_efetch", "results/discovery/candidate_accessions.tsv"])

with metadata_out.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["accession", "status", "source_genbank", "selected_from"])
    writer.writerows(metadata_rows)
