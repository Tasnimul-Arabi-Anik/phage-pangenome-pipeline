import csv
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


def remote_blast_accessions(query_fasta: Path, config):
    selection = config["reference_selection"]
    blast_output = Path(snakemake.output.hits).with_suffix(".blast6.tsv")
    cmd = [
        "blastn",
        "-query",
        str(query_fasta),
        "-db",
        selection["remote_blast_db"],
        "-remote",
        "-task",
        selection["remote_blast_task"],
        "-max_target_seqs",
        str(selection["remote_max_target_seqs"]),
        "-entrez_query",
        selection["remote_entrez_query"],
        "-outfmt",
        "6 saccver stitle pident length evalue bitscore slen",
        "-out",
        str(blast_output),
    ]
    run(cmd)

    hits = []
    seen = set()
    with blast_output.open() as handle:
        for line in handle:
            accession, title, pident, length, evalue, bitscore, slen = line.rstrip("\n").split("\t", 6)
            accession_root = accession.split(".", 1)[0]
            if accession_root in seen:
                continue
            if "phage" not in title.lower() and "virus" not in title.lower():
                continue
            genome_length = int(slen)
            if genome_length < int(selection["min_genome_length"]) or genome_length > int(selection["max_genome_length"]):
                continue
            hits.append(
                {
                    "accession": accession_root,
                    "title": title,
                    "pident": pident,
                    "alignment_length": length,
                    "evalue": evalue,
                    "bitscore": bitscore,
                    "subject_length": genome_length,
                    "source": "blastn_remote",
                }
            )
            seen.add(accession_root)
            if len(hits) >= int(selection["max_candidates"]):
                break
    return hits


def local_blast_accessions(query_fasta: Path, config):
    selection = config["reference_selection"]
    blast_db = selection.get("local_blast_db", "")
    if not blast_db:
        raise ValueError("reference_selection.local_blast_db must be set for discovery_method=local_blast_db")

    blast_output = Path(snakemake.output.hits).with_suffix(".blast6.tsv")
    cmd = [
        "blastn",
        "-query",
        str(query_fasta),
        "-db",
        blast_db,
        "-task",
        selection["local_blast_task"],
        "-max_target_seqs",
        str(selection["local_blast_max_target_seqs"]),
        "-outfmt",
        "6 saccver stitle pident length evalue bitscore slen",
        "-out",
        str(blast_output),
    ]
    run(cmd)

    hits = []
    seen = set()
    with blast_output.open() as handle:
        for line in handle:
            accession, title, pident, length, evalue, bitscore, slen = line.rstrip("\n").split("\t", 6)
            accession_root = accession.split(".", 1)[0]
            if accession_root in seen:
                continue
            genome_length = int(slen)
            if genome_length < int(selection["min_genome_length"]) or genome_length > int(selection["max_genome_length"]):
                continue
            hits.append(
                {
                    "accession": accession_root,
                    "title": title,
                    "pident": pident,
                    "alignment_length": length,
                    "evalue": evalue,
                    "bitscore": bitscore,
                    "subject_length": genome_length,
                    "source": "blastn_local_db",
                }
            )
            seen.add(accession_root)
            if len(hits) >= int(selection["max_candidates"]):
                break
    return hits


hits_out = Path(snakemake.output.hits)
candidates_out = Path(snakemake.output.candidates)
for path in [hits_out.parent, candidates_out.parent]:
    path.mkdir(parents=True, exist_ok=True)

config_inputs = snakemake.config.get("inputs", {})
source_accessions_value = config_inputs.get("existing_reference_accessions", "")
source_accessions = Path(source_accessions_value) if source_accessions_value else None
selection = snakemake.config["reference_selection"]
if selection["discovery_method"] == "configured" and source_accessions and source_accessions.exists():
    hits = [
        {
            "accession": accession.split(".", 1)[0],
            "title": "configured local cohort",
            "pident": "",
            "alignment_length": "",
            "evalue": "",
            "bitscore": "",
            "subject_length": "",
            "source": "config.inputs.existing_reference_accessions",
        }
        for accession in parse_accessions(source_accessions)
    ]
elif selection["discovery_method"] == "blastn_remote":
    hits = remote_blast_accessions(Path(snakemake.input.query), snakemake.config)
elif selection["discovery_method"] == "local_blast_db":
    hits = local_blast_accessions(Path(snakemake.input.query), snakemake.config)
else:
    raise ValueError(
        "reference_selection.discovery_method must be one of: configured, blastn_remote, local_blast_db"
    )

with hits_out.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["accession", "title", "pident", "alignment_length", "evalue", "bitscore", "subject_length", "source"])
    for hit in hits:
        writer.writerow(
            [
                hit["accession"],
                hit["title"],
                hit["pident"],
                hit["alignment_length"],
                hit["evalue"],
                hit["bitscore"],
                hit["subject_length"],
                hit["source"],
            ]
        )

with candidates_out.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["accession", "reason"])
    for hit in hits:
        if hit["source"].startswith("config"):
            reason = "configured_local_reference_set"
        elif hit["source"] == "blastn_local_db":
            reason = "blastn_local_db_candidate"
        else:
            reason = "blastn_remote_candidate"
        writer.writerow([hit["accession"], reason])
