

rule get_surf_label_from_cifti_atlas:
    input:
        cifti=lambda wildcards: config["atlas"][wildcards.atlas]["dlabel"],
    output:
        label_left="resources/atlas/atlas-{atlas}_hemi-L_parc.label.gii",
        label_right="resources/atlas/atlas-{atlas}_hemi-R_parc.label.gii",
    container:
        config["singularity"]["diffparc"]
    group:
        "grouped_{atlas}"
    shell:
        "wb_command -cifti-separate {input.cifti} COLUMN -label CORTEX_LEFT {output.label_left} -label CORTEX_RIGHT {output.label_right}"


rule resample_labels_from_fsaverage_to_fsLR:
    input:
        label="resources/atlas/atlas-{atlas}_hemi-{hemi}_parc.label.gii",
        sphere=lambda wildcards: config['template_fsavg_surf'].format(surf='sphere',hemi=wildcards.hemi),
        new_sphere=lambda wildcards: config['template_surf'].format(surf='sphere',hemi=wildcards.hemi)
    output:
        label="resources/atlas/atlas-{atlas}_space-fsLR_hemi-{hemi}_parc.label.gii",
    container:
        config["singularity"]["diffparc"]
    group:
        "grouped_subject"
    shell:
        #"wb_command -label-resample {input.label} {input.sphere} {input.new_sphere} BARYCENTRIC {output.label}"
        "cp {input.label} {output.label}"#wb_command -label-resample {input.label} {input.sphere} {input.new_sphere} BARYCENTRIC {output.label}"


rule resample_t1_surf_to_fsLR:
    input:
        surf=lambda wildcards: config["input_path"]["t1_surf"][wildcards.dataset].format(**wildcards),
        sphere=lambda wildcards: config["input_path"]["msm_sphere"][wildcards.dataset].format(**wildcards),
        new_sphere=lambda wildcards: config['template_surf'].format(surf='sphere',hemi=wildcards.hemi)
    output:
        surf=bids(
            root=root,
            datatype="anat",
            hemi="{hemi}",
            den="32k",
            suffix="{surf}.surf.gii",
            **config["subj_wildcards"],
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "grouped_subject"
    shell:
        "wb_command -surface-resample {input.surf} {input.sphere} {input.new_sphere} BARYCENTRIC {output.surf}"



rule map_atlas_to_dwi:
    input:
        vol_ref=config["input_path"]["dwi_mask"],
        label="resources/atlas/atlas-{atlas}_space-fsLR_hemi-{hemi}_parc.label.gii",
        mid_surf=bids(
            root=root,
            datatype="anat",
            hemi="{hemi}",
            den="32k",
            suffix="midthickness.surf.gii",
            **config["subj_wildcards"],
        ),
        white_surf=bids(
            root=root,
            datatype="anat",
            hemi="{hemi}",
            den="32k",
            suffix="white.surf.gii",
            **config["subj_wildcards"],
        ),
        pial_surf=bids(
            root=root,
            datatype="anat",
            hemi="{hemi}",
            den="32k",
            suffix="pial.surf.gii",
            **config["subj_wildcards"],
        ),
    output:
        vol=bids(
            root=root,
            datatype="dwi",
            hemi="{hemi}",
            atlas="{atlas}",
            suffix="dseg.nii.gz",
            **config["subj_wildcards"],
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "grouped_subject"
    shell:
        "wb_command -label-to-volume-mapping {input.label} {input.mid_surf} {input.vol_ref} {output.vol} -ribbon-constrained {input.white_surf} {input.pial_surf}"


rule resample_subcort_dseg_to_dwi:
    input:
        dseg=lambda wildcards: config["atlas"][wildcards.atlas]["label_volume"],
        vol_ref=config["input_path"]["dwi_mask"],
        warp=lambda wildcards: config["input_path"]["warp_mni_t1w"][
            wildcards.dataset
        ],
    output:
        dseg=bids(
            root=root,
            datatype="dwi",
            atlas="{atlas}",
            desc="subcort",
            space="dwi",
            suffix="dseg.nii.gz",
            **config["subj_wildcards"],
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "grouped_subject"
    shell:
        "antsApplyTransforms -i {input.dseg} -r {input.vol_ref} -t {input.warp} -o {output.dseg} -n NearestNeighbor -v"



rule merge_lr_dseg:
    input:
        left=bids(
            root=root,
            datatype="dwi",
            hemi="L",
            atlas="{atlas}",
            suffix="dseg.nii.gz",
            **config["subj_wildcards"],
        ),
        right=bids(
            root=root,
            datatype="dwi",
            hemi="R",
            atlas="{atlas}",
            suffix="dseg.nii.gz",
            **config["subj_wildcards"],
        ),
        dseg=bids(
            root=root,
            datatype="dwi",
            atlas="{atlas}",
            desc="subcort",
            space="dwi",
            suffix="dseg.nii.gz",
            **config["subj_wildcards"],
        ),
    output:
        merged=bids(
            root=root,
            datatype="dwi",
            atlas="{atlas}",
            suffix="dseg.nii.gz",
            **config["subj_wildcards"],
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "grouped_subject"
    shell:
        "c3d {input} -accum -max -endaccum -o {output.merged}"


rule get_label_txt_from_cifti:
    input:
        cifti=lambda wildcards: config["atlas"][wildcards.atlas]["dlabel"],
    output:
        label_txt="resources/atlas/atlas-{atlas}_desc-cifti_labels.txt",
    container:
        config["singularity"]["diffparc"]
    group:
        "grouped_{atlas}"
    shell:
        "wb_command -cifti-label-export-table {input} 1 {output}"


rule lut_cifti_to_bids:
    input:
        label_txt=rules.get_label_txt_from_cifti.output.label_txt,
    output:
        label_tsv=temp("resources/atlas/atlas-{atlas}_desc-nometadata_dseg.tsv"),
    group:
        "grouped_{atlas}"
    script:
        "../scripts/lut_cifti_to_bids.py"


rule lut_bids_to_itksnap:
    input:
        tsv=rules.lut_cifti_to_bids.output.label_tsv,
    output:
        lut="resources/atlas/atlas-{atlas}_desc-itksnap_labels.txt",
    group:
        "grouped_{atlas}"
    script:
        "../scripts/lut_bids_to_itksnap.py"


rule parcellate_centroids:
    input:
        dlabel=lambda wildcards: config["atlas"][wildcards.atlas]["dlabel"],
        surfs=lambda wildcards: expand(
            config["template_surf"], surf="midthickness", hemi=["L", "R"]
        ),
    params:
        method="MEDIAN",  #medoid vertex -- change to MEAN if want centroid
    output:
        coords_pscalar="resources/atlas/atlas-{atlas}_coords.pscalar.nii",
    shadow:
        "minimal"
    container:
        config["singularity"]["diffparc"]
    group:
        "grouped_{atlas}"
    shell:
        "wb_command -surface-coordinates-to-metric {input.surfs[0]} left_coords.shape.gii && "
        "wb_command -surface-coordinates-to-metric {input.surfs[1]} right_coords.shape.gii && "
        "wb_command -cifti-create-dense-scalar coords.dscalar.nii -left-metric left_coords.shape.gii -right-metric right_coords.shape.gii && "
        "wb_command -cifti-parcellate coords.dscalar.nii {input.dlabel} COLUMN {output.coords_pscalar} -method {params.method}"


rule add_metadata_to_dseg_tsv:
    input:
        coords_pscalar="resources/atlas/atlas-{atlas}_coords.pscalar.nii",
        label_tsv=rules.lut_cifti_to_bids.output.label_tsv,
    params:
        network_pattern=lambda wildcards: config["atlas"][wildcards.atlas].get(
            "network_pattern", None
        ),
    output:
        label_tsv="resources/atlas/atlas-{atlas}_dseg.tsv",
    group:
        "grouped_{atlas}"
    script:
        "../scripts/add_metadata_to_dseg_tsv.py"
