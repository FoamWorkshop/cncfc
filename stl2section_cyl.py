#!/usr/bin/python3

__author__ = 'FoamWorkshop'

import numpy as np
import stl
import cncfclib as cf
import pickle
import argparse

def prog(args):
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
    i_fname, o_fname, pedestal_params, split_list, Num_W = args
    mesh = stl.mesh.Mesh.from_file(i_fname)
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
        pos = cf.add_pedestal(pos, pedestal_params)
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
        if not o_fname:
            o_fname = i_fname

        fname='{}_{}.pickle'.format(o_fname, sect_num)
        with open(fname, 'wb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            pickle.dump(res_dict, f, pickle.HIGHEST_PROTOCOL)

        print(fname, ' saved')

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='This program converts stl files to paths for CNCFoamCutter. slender models eg. fusselage modeled in OpenVSP. The rotation axis must be along Z.')
    parser.add_argument('-i',  '--input', type=str, required=True, help='input stl file name')
    parser.add_argument('-o', '--output', type=str, help='output pickle database')
    # parser.add_argument('-sh', '--subset_header', action='store_true')
    # parser.add_argument('-gh', '--global_header', action='store_true')
    # parser.add_argument('-sw', '--swing_cut', action='store_true')
    # parser.add_argument('-cm', '--center_model', action='store_true')
    parser.add_argument('-p', '--pedestal_params', type=float, nargs='+', default=[35, 5, 5], help='pedestal dimensions')
    parser.add_argument('-s', '--sections', type = float, nargs='+', default=[0, 1000], help='location of the split sections')
    parser.add_argument('-w', '--num_w'   , type = int, default=17, help='number of longitudinal splits')
    #parser.add_argument('-symm', '--symmetry', action='store_true')

    args = parser.parse_args()

    i_fname = args.input
    o_fname = args.output
    #[radius, heigth1, heigth2]
    pedestal_params = args.pedestal_params
    split_list = args.sections
    Num_W = args.num_w

    args=(i_fname, o_fname, pedestal_params, split_list, Num_W)

    prog(args)
