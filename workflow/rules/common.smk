def get_dwi_targets():
    targets = []
    for dataset in config["datasets"]:
        targets.extend(
            expand(
                bids(
                    root=root,
                    datatype="dwi",
                    den="32k",
                    atlas="{atlas}",
                    suffix="{struc}.pconn.{plottype}.png",
                    **config["subj_wildcards"],
                ),
                subject=subjects[dataset],
                dataset=dataset,
                struc=["struc","strucFA"],
                plottype=["matrix", "chord"],
                atlas=config["atlases"],
            )
        )
        targets.extend(
            expand(
                bids(
                    root=root,
                    datatype="dwi",
                    den="32k",
                    atlas="{atlas}",
                    suffix="net{struc}.pconn.nii",
                    **config["subj_wildcards"],
                ),
                subject=subjects[dataset],
                dataset=dataset,
                struc=["struc","strucFA"],
                atlas=config["atlases"],
            )
        )
    return targets


def get_func_targets():
    targets = []
    for dataset in config["datasets"]:
        targets.extend(
            expand(
                bids(
                    root=root,
                    datatype="func",
                    desc="preproc",
                    den="32k",
                    task="{task}",
                    denoise="{denoise}",
                    fwhm="{fwhm}",
                    atlas="{atlas}",
                    suffix="bold.pconn.{plottype}.png",
                    **config["subj_wildcards"],
                ),
                subject=subjects[dataset],
                dataset=dataset,
                task=config["func"]["task"],
                denoise=config["func"]["denoise"].keys(),
                fwhm=config["func"]["fwhm"],
                atlas=config["atlases"],
                plottype=["matrix", "chord"],
            )
        )
        targets.extend(
            expand(
                bids(
                    root=root,
                    datatype="func",
                    desc="preproc",
                    den="32k",
                    task="{task}",
                    denoise="{denoise}",
                    fwhm="{fwhm}",
                    atlas="{atlas}",
                    suffix="netbold.pconn.nii",
                    **config["subj_wildcards"],
                ),
                subject=subjects[dataset],
                dataset=dataset,
                task=config["func"]["task"],
                denoise=config["func"]["denoise"].keys(),
                fwhm=config["func"]["fwhm"],
                atlas=config["atlases"],
            )
        )
    return targets


def get_sfc_targets():
    targets = []
    for dataset in config["datasets"]:
        targets.extend(
            expand(
                bids(
                    root=root,
                    datatype="func",
                    desc="preproc",
                    den="32k",
                    task="{task}",
                    denoise="{denoise}",
                    fwhm="{fwhm}",
                    atlas="{atlas}",
                    suffix="{suffix}.nii",
                    **config["subj_wildcards"],
                ),
                subject=subjects[dataset],
                dataset=dataset,
                suffix=['sfc.pscalar','netsfc.pconn'],
                task=config["func"]["task"],
                denoise=config["func"]["denoise"].keys(),
                fwhm=config["func"]["fwhm"],
                atlas=config["atlases"],
            )
        )

    return targets

def get_ttest_targets():
    targets = []
    for groups_dict in config['ttests']:
        
        targets.extend(expand(
            'stats/group1-{group1}_group2-{group2}_coltype-{coltype}_ttestplots',
                group1=groups_dict['group1'],
                group2=groups_dict['group2'],
                coltype=config['coltypes'],
                ))
    return targets
