rule qc_and_dereplicate_references:
    input:
        gb=f"results/{config['project_name']}/references/references.gb",
        metadata=f"results/{config['project_name']}/references/reference_metadata.tsv"
    output:
        close_refs=f"results/{config['project_name']}/cohorts/close_refs.txt",
        expanded_refs=f"results/{config['project_name']}/cohorts/expanded_refs.txt",
        dropped=f"results/{config['project_name']}/cohorts/dropped_references.tsv",
        cohort_meta=f"results/{config['project_name']}/cohorts/cohort_metadata.tsv"
    message:
        "Filter unusable GenBank records and derive close/expanded analysis cohorts."
    script:
        "../scripts/qc_references.py"
