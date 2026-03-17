rule plot_pangenome:
    input:
        presence_absence=f"results/{config['project_name']}/orthology/presence_absence.tsv",
        summary=f"results/{config['project_name']}/orthology/summary.tsv",
        genome_metadata=f"results/{config['project_name']}/orthology/genome_metadata.tsv"
    output:
        png=f"results/{config['project_name']}/plots/pangenome_presence_absence_heatmap.png",
        tiff=f"results/{config['project_name']}/plots/pangenome_presence_absence_heatmap.tiff"
    message:
        "Generate the orthogroup presence/absence heatmap."
    shell:
        r"""
        mkdir -p results/{config[project_name]}/plots
        python plot_pangenome_heatmap.py \
          --presence-absence {input.presence_absence} \
          --genome-metadata {input.genome_metadata} \
          --output-prefix results/{config[project_name]}/plots/pangenome_presence_absence_heatmap \
          --title "{config[plotting][title]}"
        """
