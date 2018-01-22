import numpy as np
from stl import mesh
import cncfclib as cf
import pickle
import argparse
# def main():
"""
This program converts stl files to paths for CNCFoamCutter.
slender models eg. fusselage modeled in OpenVSP.
The rotation axis must be along Z.
Parameters:
Returns:
TODO:
autodetection of the Num_W
nonblocking plots
"""
split_list = np.array([0, 120, 250])
Num_W = 17
mesh = mesh.Mesh.from_file('fuselage.stl')
#rotate mesh since by default the rotation axis is along X
mesh.rotate([0,1,0],np.pi/2)

v_arr = np.round(np.vstack(mesh.vectors).astype(float), decimals=1)

splt0_arr = np.array(split_list)
splt1_arr = np.roll(splt0_arr,-1)

pos = cf.cartesian2cylyndrical(v_arr, Num_W)

#make splits
pos_list=[]
for splt0, splt1 in zip(splt0_arr[:-1], splt1_arr[:-1]):
    pos_idx = np.where((splt0<=pos[:,:,2]) & (splt1>pos[:,:,2]))[0]
    print(splt0, splt1)
    #pos = [r, th, z] sectionwise
    pos_list.append(pos[pos_idx])
    #add pedestal mesh

for sect_num, pos in enumerate(pos_list):
    pos = cf.add_pedestal(pos)
    profiles=np.zeros_like(pos)
    for i in np.arange(np.shape(pos)[0]):
        profiles[i] = cf.cylyndrical2cartesian(pos[i])


    strokes = np.flipud(np.rot90(profiles))
    #transform data from longeron nodes [xyz] to:
    #a_arr - rotation angle around the rotation axis
    #r_arr - length of a segment perpenticular to the rotation axis and corresponding lateral mesh edge
    #z_arr - corresponding z coordiantes
    #v_arr - direction vector of the coresponding lateral mesh edge
    a_arr, r_arr, z_arr, v_arr = cf.transform(strokes, add_pedestal_bottom=True,add_pedestal_top=True)

    #make a summary plots
    cf.plot_loft_paths(profiles)
    cf.plot_loft_paths(pos)
    cf.plot_surf(a_arr,z_arr,r_arr)

    #collect data to the dictionary longeron wise
    res_dict = {'a_arr':np.rot90(a_arr, k=-1),
                'r_arr':np.rot90(r_arr, k=-1),
                'z_arr':np.rot90(z_arr, k=-1),
                'v_arr':np.rot90(v_arr, k=-1)}

    #save result dictionary
    fname='res_dict_{}.pickle'.format(sect_num)
    with open(fname, 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(res_dict, f, pickle.HIGHEST_PROTOCOL)

    print(fname, ' saved')

# if __name__ == "__main__":
#     main()
