import pandas as pd
import nibabel as nib
import numpy as np

def get_network_vertices(parcel_axis, networks):
    """ this gets the vertex indices for the agglomerated
        network regions, from the regional parcel_axis, and mapping from
        region to network"""
    structures = parcel_axis.nvertices.keys()

    network_names = networks.unique()

    verts_list = []
    for i,net in enumerate(network_names):

        region_indices = np.where(networks == net)[0]
        
        verts={}
        for structure in structures:
            
            temp_verts = []
            
            for region in region_indices:
                if structure in parcel_axis.vertices[region]:
                    temp_verts.append(parcel_axis.vertices[region][structure])
                    
            verts[structure] = np.concatenate(temp_verts)

        verts_list.append(verts)
            

    return np.array(verts_list)



nib_conn = nib.load(snakemake.input.pconn)

conn = nib_conn.get_fdata()

#get the label tsv that has roi to network lookup
df_atlas = pd.read_csv(snakemake.input.label_tsv,sep='\t',index_col='index')
df_atlas

#get list of network names
network_names = df_atlas.networks.unique()
network_names

#make the output conn matrix (networks by networks)
conn_network = np.zeros((len(network_names),len(network_names)))

#fill it in
for i,net_i in enumerate(network_names):
    for j,net_j in enumerate(network_names):
        #make a mask for the rows and cols we want (associated with each network)
        mask_rows=np.zeros(conn.shape)
        mask_cols=np.zeros(conn.shape)
        
        mask_rows[df_atlas.networks == net_i,:] = 1
        mask_cols[:,df_atlas.networks == net_j] = 1

        #get the diagonal entries as we want to remove these from the mask
        mask_diag=np.eye(len(conn))
        
        #create the final mask as a boolean
        mask_network = mask_rows * mask_cols * (mask_diag==0)
        mask_network = mask_network.astype(bool)
        
        #mask the net_i x net_j values, and take mean
        conn_vals = conn[mask_network]
        conn_network[i,j] = conn_vals.mean()


# save as a pconn cifti
# we need to update the list of regions (to be network names now):
# and also provide an updated list of vertices for each region
parcel_axis = nib_conn.header.get_axis(0)
network_axis = parcel_axis
network_axis.name = network_names
network_axis.vertices = get_network_vertices(parcel_axis,df_atlas.networks)

nib_network = nib.Cifti2Image(conn_network,header=(network_axis,network_axis))
nib_network.to_filename(snakemake.output.pconn)
