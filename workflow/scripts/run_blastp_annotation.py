import csv
import subprocess
from pathlib import Path


def write_empty_hits(path: Path):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "qseqid",
                "sseqid",
                "pident",
                "qcovs",
                "evalue",
                "bitscore",
                "subject_title",
                "blastp_product",
                "blastp_module",
            ]
        )


def write_summary(path: Path, pairs):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key, value in pairs:
            writer.writerow([key, value])


def is_query_like_subject(subject_id: str) -> bool:
    return subject_id.startswith("query|") or subject_id.startswith("query")


def load_metadata(path: Path):
    if not path.exists():
        return {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    metadata = {}
    for row in rows:
        for key in ("seq_id", "protein_id", "sseqid"):
            value = (row.get(key) or "").strip()
            if value:
                metadata[value] = row
    return metadata


query_faa = Path(snakemake.input.query_proteins)
hits_out = Path(snakemake.output.hits)
summary_out = Path(snakemake.output.summary)
hits_out.parent.mkdir(parents=True, exist_ok=True)

annotation_cfg = snakemake.config.get("annotation", {})
enabled = bool(annotation_cfg.get("blastp", False))
db_fasta_value = annotation_cfg.get("protein_db_fasta", "")
db_fasta = Path(db_fasta_value) if db_fasta_value else None
metadata_value = annotation_cfg.get("protein_db_metadata", "")
metadata_path = Path(metadata_value) if metadata_value else None

if not enabled or not db_fasta_value:
    write_empty_hits(hits_out)
    write_summary(
        summary_out,
        [
            ("enabled", str(enabled).lower()),
            ("database_configured", "false"),
            ("database_fasta", db_fasta_value),
            ("hits_retained", 0),
            ("status", "disabled_or_unconfigured"),
        ],
    )
    raise SystemExit(0)

if not db_fasta.exists():
    raise FileNotFoundError(f"Configured annotation.protein_db_fasta does not exist: {db_fasta}")

metadata = load_metadata(metadata_path) if metadata_path and metadata_path.exists() else {}

tmp_dir = hits_out.parent / "_tmp_blastp_annotation"
tmp_dir.mkdir(parents=True, exist_ok=True)
db_prefix = tmp_dir / "annotation_db"
raw_hits = tmp_dir / "query_vs_annotation_db.blast6.tsv"

makeblastdb_cmd = [
    "makeblastdb",
    "-in",
    str(db_fasta),
    "-dbtype",
    "prot",
    "-out",
    str(db_prefix),
]
subprocess.run(makeblastdb_cmd, check=True)

blastp_cmd = [
    "blastp",
    "-query",
    str(query_faa),
    "-db",
    str(db_prefix),
    "-outfmt",
    "6 qseqid sseqid pident length qlen slen qcovs evalue bitscore stitle",
    "-max_target_seqs",
    str(annotation_cfg.get("max_target_seqs", 5)),
    "-evalue",
    str(annotation_cfg.get("max_evalue", "1e-5")),
    "-out",
    str(raw_hits),
]
subprocess.run(blastp_cmd, check=True)

min_identity = float(annotation_cfg.get("min_identity", 30))
min_qcov = float(annotation_cfg.get("min_query_coverage", 50))

best_hits = {}
raw_count = 0
retained_count = 0
with raw_hits.open() as handle:
    reader = csv.reader(handle, delimiter="\t")
    for row in reader:
        if len(row) < 10:
            continue
        raw_count += 1
        qseqid, sseqid, pident, _length, _qlen, _slen, qcovs, evalue, bitscore, stitle = row
        if is_query_like_subject(sseqid):
            continue
        pident_f = float(pident)
        qcovs_f = float(qcovs)
        evalue_f = float(evalue)
        bitscore_f = float(bitscore)
        if pident_f < min_identity or qcovs_f < min_qcov:
            continue
        retained_count += 1
        meta = metadata.get(sseqid) or metadata.get(sseqid.split("|")[-1], {})
        meta_product = (meta.get("product") or "").strip()
        meta_module = (meta.get("module") or "").strip()
        product = meta_product or stitle.strip() or sseqid
        current = best_hits.get(qseqid)
        candidate = {
            "qseqid": qseqid,
            "sseqid": sseqid,
            "pident": pident_f,
            "qcovs": qcovs_f,
            "evalue": evalue_f,
            "bitscore": bitscore_f,
            "subject_title": meta_product or stitle.strip(),
            "blastp_product": product,
            "blastp_module": meta_module,
        }
        if current is None or (
            candidate["bitscore"],
            -candidate["evalue"],
            candidate["pident"],
        ) > (
            current["bitscore"],
            -current["evalue"],
            current["pident"],
        ):
            best_hits[qseqid] = candidate

with hits_out.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "qseqid",
            "sseqid",
            "pident",
            "qcovs",
            "evalue",
            "bitscore",
            "subject_title",
            "blastp_product",
            "blastp_module",
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
        ("database_fasta", db_fasta),
        ("database_metadata", metadata_path or ""),
        ("raw_hits", raw_count),
        ("hits_passing_thresholds", retained_count),
        ("best_hits_written", len(best_hits)),
        ("min_identity", min_identity),
        ("min_query_coverage", min_qcov),
        ("max_evalue", annotation_cfg.get("max_evalue", "1e-5")),
        ("status", "completed"),
    ],
)
