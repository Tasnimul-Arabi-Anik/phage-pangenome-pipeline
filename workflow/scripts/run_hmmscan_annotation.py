import csv
import subprocess
from pathlib import Path


def write_empty_hits(path: Path):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "qseqid",
                "pfam_name",
                "pfam_accession",
                "pfam_evalue",
                "pfam_score",
                "pfam_description",
                "pfam_product",
            ]
        )


def write_summary(path: Path, pairs):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key, value in pairs:
            writer.writerow([key, value])


query_faa = Path(snakemake.input.query_proteins)
hits_out = Path(snakemake.output.hits)
summary_out = Path(snakemake.output.summary)
hits_out.parent.mkdir(parents=True, exist_ok=True)

annotation_cfg = snakemake.config.get("annotation", {})
enabled = bool(annotation_cfg.get("pfam_hmmscan", False))
pfam_db_value = annotation_cfg.get("pfam_db", "")
pfam_db = Path(pfam_db_value) if pfam_db_value else None

if not enabled or not pfam_db_value:
    write_empty_hits(hits_out)
    write_summary(
        summary_out,
        [
            ("enabled", str(enabled).lower()),
            ("database_configured", "false"),
            ("pfam_db", pfam_db_value),
            ("best_hits_written", 0),
            ("status", "disabled_or_unconfigured"),
        ],
    )
    raise SystemExit(0)

if not pfam_db.exists():
    raise FileNotFoundError(f"Configured annotation.pfam_db does not exist: {pfam_db}")

tblout_path = hits_out.parent / "_tmp_query_pfam.tblout"
hmmscan_cmd = [
    "hmmscan",
    "--tblout",
    str(tblout_path),
    str(pfam_db),
    str(query_faa),
]
subprocess.run(hmmscan_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

max_evalue = float(annotation_cfg.get("pfam_max_evalue", 1e-3))
best_hits = {}
raw_hits = 0

with tblout_path.open() as handle:
    for raw_line in handle:
        if not raw_line.strip() or raw_line.startswith("#"):
            continue
        raw_hits += 1
        parts = raw_line.strip().split(maxsplit=18)
        if len(parts) < 18:
            continue
        target_name = parts[0]
        target_accession = parts[1]
        query_name = parts[2]
        evalue = float(parts[4])
        score = float(parts[5])
        description = parts[18] if len(parts) > 18 else target_name
        if evalue > max_evalue:
            continue
        candidate = {
            "qseqid": query_name,
            "pfam_name": target_name,
            "pfam_accession": target_accession,
            "pfam_evalue": evalue,
            "pfam_score": score,
            "pfam_description": description,
            "pfam_product": description or target_name,
        }
        current = best_hits.get(query_name)
        if current is None or (candidate["pfam_score"], -candidate["pfam_evalue"]) > (
            current["pfam_score"],
            -current["pfam_evalue"],
        ):
            best_hits[query_name] = candidate

with hits_out.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "qseqid",
            "pfam_name",
            "pfam_accession",
            "pfam_evalue",
            "pfam_score",
            "pfam_description",
            "pfam_product",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    for qseqid in sorted(best_hits):
        writer.writerow(best_hits[qseqid])

write_summary(
    summary_out,
    [
        ("enabled", "true"),
        ("database_configured", "true"),
        ("pfam_db", pfam_db),
        ("raw_hits", raw_hits),
        ("best_hits_written", len(best_hits)),
        ("max_evalue", max_evalue),
        ("status", "completed"),
    ],
)
