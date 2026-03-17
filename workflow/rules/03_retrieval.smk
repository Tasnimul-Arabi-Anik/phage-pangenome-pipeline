rule retrieve_references:
    input:
        candidates=f"results/{config['project_name']}/discovery/candidate_accessions.tsv"
    output:
        gb=f"results/{config['project_name']}/references/references.gb",
        fna=f"results/{config['project_name']}/references/references.fna",
        metadata=f"results/{config['project_name']}/references/reference_metadata.tsv"
    message:
        "Stage the configured reference GenBank file into the workflow outputs."
    script:
        "../scripts/retrieve_references.py"
