rule validate_input:
    input:
        query=lambda wildcards: config["query_fasta"]
    output:
        manifest=f"results/{config['project_name']}/input/input_manifest.tsv"
    message:
        "Validate query FASTA and initialize the project manifest."
    script:
        "../scripts/validate_input.py"
