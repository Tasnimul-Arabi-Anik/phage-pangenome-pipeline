rule classify_query_genes:
    input:
        presence_absence=f"results/{config['project_name']}/orthology/presence_absence.tsv",
        manifest=f"results/{config['project_name']}/proteins/protein_manifest.tsv"
    output:
        query_table=f"results/{config['project_name']}/interpretation/query_gene_orthogroup_classification.tsv",
        query_singletons=f"results/{config['project_name']}/interpretation/query_singletons.tsv"
    message:
        "Classify query genes into core, accessory, and singleton orthogroups."
    script:
        "../scripts/classify_query_genes.py"
