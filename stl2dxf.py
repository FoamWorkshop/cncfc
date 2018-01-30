#!/usr/bin/python3

__author__ = 'FoamWorkshop'

'''The program automaticaly extracts cutting path from a dxf file.
The cutting path is split into:
1. IO_path - in/out path begining with single knot
2. ct_path - closed loop, begin and end in master knot position
the output of the program is a set of files with an ordered list of knots'''

import numpy as np
from stl import mesh
import cncfclib as fc
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from dxfwrite import DXFEngine as dxf
import argparse

parser = argparse.ArgumentParser(description='test')
parser.add_argument('-i',  '--input', default=[], type=str, help='input stl file')
parser.add_argument('-o', '--output', default=[], type=str, help='output dxf file')
parser.add_argument('-s', '--slices', default=[10,1,1], nargs='+', type=int, help='number of slices per axis')

args = parser.parse_args()
fname_stl = args.input
fname_dxf = args.output
slices_list = args.slices

if ~len(fname_dxf):
    fname_dxf=fname_stl.split('.')[0]+'.dxf'

mesh = mesh.Mesh.from_file(fname_stl)

dim_max = np.round(mesh.max_,3)
dim_min = np.round(mesh.min_,3)

sections_list = [([1,0,0], dim_min[0], dim_max[0], slices_list[0]),
                 ([0,1,0], dim_min[1], dim_max[1], slices_list[1]),
                 ([0,0,1], dim_min[2], dim_max[2], slices_list[2])]

section_plane_list=[]
for i, (norm, n_min, n_max, n_sect) in enumerate(sections_list):
    print('slicing plane: ',i)
    # n_sect=10
    cp_n0_arr = np.tile(norm,(n_sect,1))

    if n_sect==1:
        cp_D0_arr = np.tile(norm,(n_sect,1)) * np.array([0.5*(n_min + n_max)])
    else:
        cp_D0_arr = np.tile(norm,(n_sect,1)) * np.linspace(n_min, n_max, n_sect)[:,np.newaxis]

    section_list=[]
    #
    for i, (n0, D0) in enumerate(zip(cp_n0_arr, cp_D0_arr)):
        intersect_list = []
        for tri_list in mesh.vectors:
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
                    block_x.add(dxf.line((z0+x0, y0), (z1+x1, y1),thickness=0, color=j))
                    # block_x.add(dxf.line((y0+x0, z0), (y1+x1, z1),thickness=0, color=j))
                    block_x_ref = dxf.insert(blockname_x, insert=(0,0), layer='~slices_x')
                    drawing.add(block_x_ref)

                elif i==1:
                    block_y.add(dxf.line((x0, z0+y0), (x1, z1+y0),thickness=0, color=j))
                    # block_y.add(dxf.line((x0, z0), (x1, z1),thickness=0, color=j))
                    block_y_ref = dxf.insert(blockname_y, insert=(0,0), layer='~slices_y')
                    drawing.add(block_y_ref)

                elif i==2:
                    block_z.add(dxf.line((x0, y0), (x1, y1),thickness=0, color=j))
                    block_z_ref = dxf.insert(blockname_z, insert=(0,0), layer='~slices_z')
                    drawing.add(block_z_ref)

        #
drawing.save()


# fig = plt.figure()
# ax = fig.gca(projection='3d')
#
# for section_list in section_plane_list:
#     for i, section in enumerate(section_list):
#         p_arr =  np.array(section)
#         for row in p_arr:
#             if row.shape[0]==2:
#                 # print(row)
#                 x = row[:,0]
#                 y = row[:,1]
#                 z = row[:,2]
#                 ax.plot(x, y, z,'s-')
# plt.show()
print('Thank you')
