import csv
import subprocess
from collections import Counter, defaultdict
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


def load_manifest(path: Path):
    rows = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows[row["seq_id"]] = row
    return rows


def load_proteins(fasta_path: Path, manifest_path: Path):
    manifest = load_manifest(manifest_path)
    proteins = []
    genomes = []
    seen = set()
    for header, sequence in parse_fasta(fasta_path):
        row = manifest[header.split()[0]]
        proteins.append(
            Protein(
                seq_id=row["seq_id"],
                genome_id=row["genome_id"],
                genome_name=row["genome_name"],
                protein_id=row["protein_id"],
                locus_tag=row["locus_tag"],
                product=row["product"],
                module=row["module"],
                sequence=sequence.rstrip("*"),
                source=row["source"],
            )
        )
        if row["genome_id"] not in seen:
            genomes.append((row["genome_id"], row["genome_name"]))
            seen.add(row["genome_id"])
    return proteins, genomes


def run(cmd):
    subprocess.run(cmd, check=True)


def blast_all_vs_all(combined_fasta: Path, blast_db_prefix: Path, blast_out: Path):
    run(["makeblastdb", "-in", str(combined_fasta), "-dbtype", "prot", "-out", str(blast_db_prefix)])
    run(
        [
            "blastp",
            "-query",
            str(combined_fasta),
            "-db",
            str(blast_db_prefix),
            "-outfmt",
            "6 qseqid sseqid pident length qlen slen evalue bitscore",
            "-seg",
            "no",
            "-out",
            str(blast_out),
        ]
    )


def filtered_hits(path: Path, min_identity: float, min_query_cov: float, min_subject_cov: float, max_evalue: float):
    hits = []
    with path.open() as handle:
        for line in handle:
            qseqid, sseqid, pident, length, qlen, slen, evalue, bitscore = line.rstrip("\n").split("\t")
            if qseqid == sseqid:
                continue
            pident = float(pident)
            length = int(length)
            qlen = int(qlen)
            slen = int(slen)
            evalue = float(evalue)
            bitscore = float(bitscore)
            qcov = 100.0 * length / qlen if qlen else 0.0
            scov = 100.0 * length / slen if slen else 0.0
            if pident < min_identity or qcov < min_query_cov or scov < min_subject_cov or evalue > max_evalue:
                continue
            hits.append(
                {
                    "qseqid": qseqid,
                    "sseqid": sseqid,
                    "pident": pident,
                    "qcov": qcov,
                    "scov": scov,
                    "bitscore": bitscore,
                }
            )
    return hits


def choose_best_hits(hits, proteins_by_id):
    best = {}
    for hit in hits:
        qseqid = hit["qseqid"]
        sseqid = hit["sseqid"]
        if proteins_by_id[qseqid].genome_id == proteins_by_id[sseqid].genome_id:
            continue
        key = (qseqid, proteins_by_id[sseqid].genome_id)
        score = (hit["bitscore"], hit["pident"], min(hit["qcov"], hit["scov"]))
        if key not in best or score > best[key][0]:
            best[key] = (score, sseqid)
    return {key: value[1] for key, value in best.items()}


def build_components(nodes, edges):
    parent = {node: node for node in nodes}

    def find(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(left, right):
        root_left = find(left)
        root_right = find(right)
        if root_left != root_right:
            parent[root_right] = root_left

    for left, right in edges:
        union(left, right)

    groups = defaultdict(list)
    for node in nodes:
        groups[find(node)].append(node)
    return list(groups.values())


def make_clusters(proteins, hits):
    proteins_by_id = {protein.seq_id: protein for protein in proteins}
    best = choose_best_hits(hits, proteins_by_id)
    edges = set()
    for (qseqid, subject_genome), sseqid in best.items():
        reciprocal = best.get((sseqid, proteins_by_id[qseqid].genome_id))
        if reciprocal == qseqid:
            edges.add(tuple(sorted((qseqid, sseqid))))
    return build_components(proteins_by_id.keys(), edges), proteins_by_id, edges


def summarize_product(values):
    cleaned = [value for value in values if value]
    informative = [value for value in cleaned if "hypothetical" not in value.lower()]
    if informative:
        return Counter(informative).most_common(1)[0][0]
    if cleaned:
        return Counter(cleaned).most_common(1)[0][0]
    return "hypothetical protein"


def classify_group(genome_count, total_genomes):
    if genome_count == total_genomes:
        return "core"
    if genome_count == 1:
        return "singleton"
    return "accessory"


def write_outputs(proteins, genomes, components, edges, output_dir: Path, summary_path: Path, config):
    proteins_by_id = {protein.seq_id: protein for protein in proteins}
    genome_order = [genome_id for genome_id, _ in genomes]

    orthogroups = []
    for idx, members in enumerate(sorted(components, key=lambda group: (-len(group), sorted(group)[0])), start=1):
        member_proteins = [proteins_by_id[member] for member in sorted(members)]
        genome_presence = Counter(protein.genome_id for protein in member_proteins)
        orthogroups.append(
            {
                "orthogroup": f"OG{idx:04d}",
                "members": member_proteins,
                "genome_presence": genome_presence,
                "category": classify_group(len(genome_presence), len(genome_order)),
                "consensus_product": summarize_product([protein.product for protein in member_proteins]),
                "query_module": summarize_product(
                    [protein.module for protein in member_proteins if protein.genome_id == "query" and protein.module]
                ),
            }
        )

    with Path(snakemake.output.orthogroups).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "orthogroup",
                "category",
                "consensus_product",
                "query_module",
                "n_genomes",
                "n_proteins",
                "members",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for orthogroup in orthogroups:
            writer.writerow(
                {
                    "orthogroup": orthogroup["orthogroup"],
                    "category": orthogroup["category"],
                    "consensus_product": orthogroup["consensus_product"],
                    "query_module": orthogroup["query_module"],
                    "n_genomes": len(orthogroup["genome_presence"]),
                    "n_proteins": len(orthogroup["members"]),
                    "members": ",".join(protein.seq_id for protein in orthogroup["members"]),
                }
            )

    with Path(snakemake.output.presence_absence).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["orthogroup", "category", "consensus_product", "query_module", "n_genomes"] + genome_order,
            delimiter="\t",
        )
        writer.writeheader()
        for orthogroup in orthogroups:
            row = {
                "orthogroup": orthogroup["orthogroup"],
                "category": orthogroup["category"],
                "consensus_product": orthogroup["consensus_product"],
                "query_module": orthogroup["query_module"],
                "n_genomes": len(orthogroup["genome_presence"]),
            }
            for genome_id in genome_order:
                members = [protein.protein_id for protein in orthogroup["members"] if protein.genome_id == genome_id]
                row[genome_id] = ",".join(members)
            writer.writerow(row)

    genome_metadata_path = Path(snakemake.output.genome_metadata)
    with genome_metadata_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["genome_id", "genome_name", "protein_count"])
        for genome_id, genome_name in genomes:
            writer.writerow([genome_id, genome_name, sum(1 for protein in proteins if protein.genome_id == genome_id)])

    category_counts = Counter(orthogroup["category"] for orthogroup in orthogroups)
    query_total = sum(1 for protein in proteins if protein.genome_id == "query")
    query_core = sum(1 for orthogroup in orthogroups if orthogroup["category"] == "core" and orthogroup["genome_presence"].get("query"))
    query_accessory = sum(1 for orthogroup in orthogroups if orthogroup["category"] == "accessory" and orthogroup["genome_presence"].get("query"))
    query_singleton = sum(1 for orthogroup in orthogroups if orthogroup["category"] == "singleton" and orthogroup["genome_presence"].get("query"))

    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["project_name", config["project_name"]])
        writer.writerow(["mode", config["mode"]])
        writer.writerow(["genomes", len(genome_order)])
        writer.writerow(["reference_genomes", len(genome_order) - 1])
        writer.writerow(["proteins_total", len(proteins)])
        writer.writerow(["orthogroups_total", len(orthogroups)])
        writer.writerow(["core_orthogroups", category_counts["core"]])
        writer.writerow(["accessory_orthogroups", category_counts["accessory"]])
        writer.writerow(["singleton_orthogroups", category_counts["singleton"]])
        writer.writerow(["query_proteins_total", query_total])
        writer.writerow(["query_core_orthogroups", query_core])
        writer.writerow(["query_accessory_orthogroups", query_accessory])
        writer.writerow(["query_singleton_orthogroups", query_singleton])
        writer.writerow(["rbh_edges", len(edges)])
        writer.writerow(["min_identity", config["orthology"]["min_identity"]])
        writer.writerow(["min_query_coverage", config["orthology"]["min_query_coverage"]])
        writer.writerow(["min_subject_coverage", config["orthology"]["min_subject_coverage"]])
        writer.writerow(["max_evalue", config["orthology"]["max_evalue"]])


combined_fasta = Path(snakemake.input.proteins)
manifest_path = Path(snakemake.input.manifest)
orthology_dir = Path(snakemake.output.orthogroups).parent
orthology_dir.mkdir(parents=True, exist_ok=True)
tmp_dir = orthology_dir / "_tmp"
tmp_dir.mkdir(parents=True, exist_ok=True)
blast_db_prefix = tmp_dir / "pangenome_db"
blast_out = tmp_dir / "all_vs_all_blastp.tsv"

proteins, genomes = load_proteins(combined_fasta, manifest_path)
blast_all_vs_all(combined_fasta, blast_db_prefix, blast_out)
hits = filtered_hits(
    blast_out,
    float(snakemake.config["orthology"]["min_identity"]),
    float(snakemake.config["orthology"]["min_query_coverage"]),
    float(snakemake.config["orthology"]["min_subject_coverage"]),
    float(snakemake.config["orthology"]["max_evalue"]),
)
components, proteins_by_id, edges = make_clusters(proteins, hits)
write_outputs(proteins, genomes, components, edges, orthology_dir, Path(snakemake.output.summary), snakemake.config)
