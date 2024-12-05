import numpy as np

def get_network_vertices_voxels(parcel_axis, networks):
    """ this gets the vertex indices for the agglomerated
        network regions, from the regional parcel_axis, and mapping from
        region to network"""
    structures = parcel_axis.nvertices.keys()

    network_names = networks.unique()

    verts_list = []
    voxels_list = []
    for i,net in enumerate(network_names):
        region_indices = np.where(networks == net)[0]
        
        verts={}
        voxels={}
        for structure in structures:

            temp_verts = []
            temp_voxels = []
            
            for region in region_indices:
                if structure in parcel_axis.vertices[region]:
                    temp_verts.append(parcel_axis.vertices[region][structure])
                       
            if len(temp_verts)>0:
                verts[structure] = np.concatenate(temp_verts)

        #voxels (at least how we have mapped them), don't seem to be associated with
        #any structure..
            
        for region in region_indices:        
            temp_voxels.append(parcel_axis.voxels[region])

        voxels_array = np.concatenate(temp_voxels)
        voxels_list.append(voxels_array)
        verts_list.append(verts)
            
    return (np.array(verts_list),voxels_list)



