rule import_hcp_tabular:
    input:
        ndar_txt=lambda wildcards: Path(in_root) / config["in_tabular"][wildcards.dataset],
    params:
        subjects=lambda wildcards: subjects[wildcards.dataset]
    output:
        tsv="resources/dataset-{dataset,HCP}_tabular.tsv",
    script:
        "../scripts/import_hcp_tabular.py"

rule import_lobe_tabular:
    input:
        tsv=lambda wildcards: Path(in_root) / config["in_tabular"][wildcards.dataset],
    params:
        subjects=lambda wildcards: subjects[wildcards.dataset]
    output:
        tsv="resources/dataset-{dataset,LOBE}_tabular.tsv",
    script:
        "../scripts/import_lobe_tabular.py"


