rule characterize_query:
    input:
        query=config["query_fasta"],
        manifest=f"results/{config['project_name']}/input/input_manifest.tsv"
    output:
        metrics=f"results/{config['project_name']}/query/query_metrics.tsv",
        proteins=f"results/{config['project_name']}/query/query_proteins.faa",
        cds=f"results/{config['project_name']}/query/query_cds.fna",
        gff=f"results/{config['project_name']}/query/query.gff"
    message:
        "Characterize the query genome and normalize the ORF/protein outputs."
    script:
        "../scripts/characterize_query.py"
