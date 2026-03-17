rule build_combined_proteins:
    input:
        query_proteins=f"results/{config['project_name']}/query/query_proteins.faa",
        gb=f"results/{config['project_name']}/references/references.gb",
        cohort_meta=f"results/{config['project_name']}/cohorts/cohort_metadata.tsv"
    output:
        combined=f"results/{config['project_name']}/proteins/combined_proteins.faa",
        manifest=f"results/{config['project_name']}/proteins/protein_manifest.tsv"
    message:
        "Extract translated reference CDS entries and harmonize them with the query proteins."
    script:
        "../scripts/build_combined_proteins.py"
