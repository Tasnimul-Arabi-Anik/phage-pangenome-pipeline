rule feature_specific_followup:
    input:
        query_table=f"results/{config['project_name']}/interpretation/query_gene_orthogroup_classification.tsv"
    output:
        note=f"results/{config['project_name']}/features/feature_followup_note.md"
    message:
        "Generate a query-centered feature follow-up summary."
    script:
        "../scripts/generate_feature_note.py"
