rule infer_orthogroups:
    input:
        proteins=f"results/{config['project_name']}/proteins/combined_proteins.faa",
        manifest=f"results/{config['project_name']}/proteins/protein_manifest.tsv"
    output:
        orthogroups=f"results/{config['project_name']}/orthology/orthogroups.tsv",
        presence_absence=f"results/{config['project_name']}/orthology/presence_absence.tsv",
        summary=f"results/{config['project_name']}/orthology/summary.tsv",
        genome_metadata=f"results/{config['project_name']}/orthology/genome_metadata.tsv"
    message:
        "Infer RBH orthogroups and write pangenome summary tables."
    script:
        "../scripts/infer_orthogroups.py"
