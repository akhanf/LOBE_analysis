import pandas as pd
import nibabel as nib
import numpy as np
from scipy.stats import spearmanr
from lib.utils import get_network_vertices_voxels



nib_sc = nib.load(snakemake.input.pconn_sc)
nib_fc = nib.load(snakemake.input.pconn_fc)


sc = nib_sc.get_fdata()
fc = nib_fc.get_fdata()

#get the label tsv that has roi to network lookup
df_atlas = pd.read_csv(snakemake.input.label_tsv,sep='\t',index_col='index')
df_atlas

#get list of network names
network_names = df_atlas.networks.unique()
network_names

#make the output conn matrix (networks by networks)
sfc_network = np.zeros((len(network_names),len(network_names)))

#fill it in
for i,net_i in enumerate(network_names):
    for j,net_j in enumerate(network_names):
        #make a mask for the rows and cols we want (associated with each network)
        mask_rows=np.zeros(sc.shape)
        mask_cols=np.zeros(sc.shape)
        
        mask_rows[df_atlas.networks == net_i,:] = 1
        mask_cols[:,df_atlas.networks == net_j] = 1

        #get the diagonal entries as we want to remove these from the mask
        mask_diag=np.eye(len(sc))
        
        #create the final mask as a boolean
        mask_network = mask_rows * mask_cols * (mask_diag==0)
        mask_network = mask_network.astype(bool)
        
        #mask the net_i x net_j values
        sc_vals = sc[mask_network]
        fc_vals = fc[mask_network]

    
        sfc_network[i,j] = spearmanr(sc_vals,fc_vals).statistic


# save as a pconn cifti
# we need to update the list of regions (to be network names now):
# and also provide an updated list of vertices for each region
parcel_axis = nib_sc.header.get_axis(0)
network_axis = parcel_axis
network_axis.name = network_names
(network_axis.vertices,network_axis.voxels) = get_network_vertices_voxels(parcel_axis,df_atlas.networks)

nib_network = nib.Cifti2Image(sfc_network,header=(network_axis,network_axis))

nib_network.to_filename(snakemake.output.pconn)
