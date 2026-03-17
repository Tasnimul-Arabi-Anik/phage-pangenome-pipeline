rule build_report:
    input:
        summary=f"results/{config['project_name']}/orthology/summary.tsv",
        query_table=f"results/{config['project_name']}/interpretation/query_gene_orthogroup_classification.tsv",
        feature_note=f"results/{config['project_name']}/features/feature_followup_note.md",
        heatmap=f"results/{config['project_name']}/plots/pangenome_presence_absence_heatmap.png"
    output:
        md=f"results/{config['project_name']}/report/report.md",
        docx=f"results/{config['project_name']}/report/report.docx"
    message:
        "Build a manuscript-style Markdown and DOCX report from the pipeline outputs."
    shell:
        r"""
        mkdir -p results/{config[project_name]}/report
        python workflow/scripts/build_report.py \
          --summary {input.summary} \
          --query-table {input.query_table} \
          --feature-note {input.feature_note} \
          --heatmap {input.heatmap} \
          --output-md {output.md} \
          --project-name {config[project_name]} \
          --mode {config[mode]}
        pandoc {output.md} -o {output.docx}
        """
