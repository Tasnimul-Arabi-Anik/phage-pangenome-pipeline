rule enrich_query_annotations:
    input:
        query_table=f"results/{config['project_name']}/interpretation/query_gene_orthogroup_classification.tsv"
    output:
        annotation_table=f"results/{config['project_name']}/features/query_gene_annotations.tsv",
        summary=f"results/{config['project_name']}/features/annotation_summary.tsv"
    message:
        "Enrich query annotations from orthogroup consensus labels."
    script:
        "../scripts/enrich_query_annotations.py"


rule feature_specific_followup:
    input:
        query_table=f"results/{config['project_name']}/features/query_gene_annotations.tsv"
    output:
        note=f"results/{config['project_name']}/features/feature_followup_note.md"
    message:
        "Generate a query-centered feature follow-up summary."
    script:
        "../scripts/generate_feature_note.py"
