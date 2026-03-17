rule discover_references:
    input:
        query=config["query_fasta"],
        metrics=f"results/{config['project_name']}/query/query_metrics.tsv"
    output:
        hits=f"results/{config['project_name']}/discovery/blast_hits.tsv",
        candidates=f"results/{config['project_name']}/discovery/candidate_accessions.tsv"
    message:
        "Record the configured local reference cohort for downstream staging."
    script:
        "../scripts/discover_references.py"
