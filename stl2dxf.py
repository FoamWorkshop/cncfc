#!/usr/bin/python3

__author__ = 'FoamWorkshop'
__version__ = 0.01
''' '''

import numpy as np
from stl import mesh
import cncfclib as fc
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from dxfwrite import DXFEngine as dxf
import argparse

def main(args):
    fname_stl = args.input
    fname_dxf = args.output
    slices_list = args.slices
    k_e = args.slices_end_factors
    proj_offsets = args.proj_offsets
    plot_sections = args.plot_sections

    # fname_stl, fname_dxf, slices_list, k_e, proj_offsets, plot_sections = args_t

    sx = args.x_slices
    sy = args.y_slices
    sz = args.z_slices
    print(sx, sy, sz)
    if not len(fname_dxf):
        fname_dxf=fname_stl.split('.')[0]+'.dxf'
    print(fname_stl)
    mesh_stl = mesh.Mesh.from_file(fname_stl)

    dim_max = np.round(mesh_stl.max_,3)
    dim_min = np.round(mesh_stl.min_,3)

    sections_list = [([1,0,0], dim_min[0] + k_e[0], dim_max[0] - k_e[0], slices_list[0], sx),
                     ([0,1,0], dim_min[1] + k_e[1], dim_max[1] - k_e[1], slices_list[1], sy),
                     ([0,0,1], dim_min[2] + k_e[2], dim_max[2] - k_e[2], slices_list[2], sz)]

    section_plane_list=[]
    for i, (norm, n_min, n_max, n_sect, s) in enumerate(sections_list):
        print('slicing plane: ',i)
        # n_sect=10
        cp_n0_arr = np.tile(norm,(n_sect,1))

        if n_sect==1:
            cp_D0_arr = np.tile(norm,(n_sect,1)) * np.array([0.5*(n_min + n_max)])
        else:
            #k_e jest offsetem, ktory zapobiega rzutowaniu sciany na plaszczyzne.
            #rzutowanie triangulowanej plaszczyzny skutkuje 'poszatkowaniem' rzutu i
            #dodatkowej, manualnej obrobki dxfa
            cp_D0_arr = np.tile(norm,(n_sect,1)) * np.linspace(n_min, n_max, n_sect)[:,np.newaxis]

        if len(s):
            s_cp_n0_arr = np.tile(norm,(len(s),1))
            s_cp_D0_arr = np.tile(norm,(len(s),1)) * np.array(s)[:,np.newaxis]
            cp_D0_arr = np.vstack([cp_D0_arr, s_cp_D0_arr])
            cp_D0_arr = np.unique(cp_D0_arr, axis=0)
            cp_n0_arr = np.tile(norm,(cp_D0_arr.shape[0],1))

            # print(np.vstack([cp_n0_arr, s_cp_n0_arr]))

        section_list=[]
        #
        for i, (n0, D0) in enumerate(zip(cp_n0_arr, cp_D0_arr)):
            intersect_list = []
            # print(i)
            for tri_list in mesh_stl.vectors:
                tri = np.round(np.vstack(tri_list).astype(float),6)
                intersect = fc.tri_plane_intersect(tri, D0, n0)
                if len(intersect)==2:
                    intersect_list.append(np.vstack(intersect))
            section_list.append(intersect_list)
        section_plane_list.append(section_list)

    drawing = dxf.drawing(fname_dxf)
    drawing.add_layer('~slices_x')
    drawing.add_layer('~slices_y')
    drawing.add_layer('~slices_z')

    for i, section_list in enumerate(section_plane_list):
        for j, section in enumerate(section_list):
            blockname_x = '_'.join(('slices_x',str(j)))
            blockname_y = '_'.join(('slices_y',str(j)))
            blockname_z = '_'.join(('slices_z',str(j)))
            block_x = dxf.block(name = blockname_x)
            block_y = dxf.block(name = blockname_y)
            block_z = dxf.block(name = blockname_z)
            drawing.add(block_x)
            drawing.add(block_y)
            drawing.add(block_z)

            p_arr =  np.array(section)
            for row in p_arr:
                if row.shape[0]==2:
                    x0 = row[0,0]
                    y0 = row[0,1]
                    z0 = row[0,2]

                    x1 = row[1,0]
                    y1 = row[1,1]
                    z1 = row[1,2]

                    if i==0:
                        if proj_offsets:
                            block_x.add(dxf.line((z0+x0, y0), (z1+x1, y1),thickness=0, color=j))
                        else:
                            block_x.add(dxf.line((y0, z0), (y1, z1),thickness=0, color=j))

                        block_x_ref = dxf.insert(blockname_x, insert=(0,0), layer='~slices_x')
                        drawing.add(block_x_ref)

                    elif i==1:
                        if proj_offsets:
                            block_y.add(dxf.line((x0, z0+y0), (x1, z1+y0),thickness=0, color=j))
                        else:
                            block_y.add(dxf.line((x0, z0), (x1, z1),thickness=0, color=j))

                        block_y_ref = dxf.insert(blockname_y, insert=(0,0), layer='~slices_y')
                        drawing.add(block_y_ref)

                    elif i==2:
                        if proj_offsets:
                            block_z.add(dxf.line((x0, y0), (x1, y1),thickness=0, color=j))
                        else:
                            block_z.add(dxf.line((x0, y0), (x1, y1),thickness=0, color=j))

                        block_z_ref = dxf.insert(blockname_z, insert=(0,0), layer='~slices_z')
                        drawing.add(block_z_ref)
    drawing.save()

    if plot_sections:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        for section_list in section_plane_list:
            for i, section in enumerate(section_list):
                p_arr =  np.array(section)
                for row in p_arr:
                    if row.shape[0]==2:
                        # print(row)
                        x = row[:,0]
                        y = row[:,1]
                        z = row[:,2]
                        ax.plot(x, y, z,'s-')
        plt.show()
    print('Thank you')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='This program slices dxf model using XYZ planes and saves projection of the section curves to a dxf file. Each section curve is saved as a block inserted into corresponding layer. In order to make the usage easier, block are colored.')
    parser.add_argument('-i', '--input', required = True, type=str, help='stl input file')
    parser.add_argument('-o', '--output', default=[], type=str, help='dxf output file')
    parser.add_argument('-s', '--slices', default=[10,1,1], nargs='+', type=int, help='(3 int numbers) number of slices per axis.\n 1 section -> mid plane\n2 - end section planes, >2 - linear spacing between min and max dimension')
    parser.add_argument('-k_e', '--slices_end_factors', default=[1,1,1], nargs='+', type=float, help='end sections offset (3 float numbers)')
    parser.add_argument('-p_o', '--proj_offsets', action='store_true', help='projection offsets. if activated additional coord ')
    parser.add_argument('-p_s', '--plot_sections', action='store_true', help='plot 3d sections using matplotlib')
    parser.add_argument('-sx', '--x_slices', default=[], nargs='+', type=float, help='define x slices, if -s is specified, the 2 lists are merged')
    parser.add_argument('-sy', '--y_slices', default=[], nargs='+', type=float, help='define y slices, if -s is specified, the 2 lists are merged')
    parser.add_argument('-sz', '--z_slices', default=[], nargs='+', type=float, help='define z slices, if -s is specified, the 2 lists are merged')
    args = parser.parse_args()
    # print(args)
    # fname_stl = args.input
    # fname_dxf = args.output
    # slices_list = args.slices
    # k_e = args.slices_end_factors
    # proj_offsets = args.proj_offsets
    # plot_sections = args.plot_sections


    main(args)

else:
    args=[]
    main(args)
