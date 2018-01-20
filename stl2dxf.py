import numpy as np
from stl import mesh
import cncfclib as fc
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from dxfwrite import DXFEngine as dxf


mesh = mesh.Mesh.from_file('fuselage.stl')

dim_max = np.round(mesh.max_,3)
dim_min = np.round(mesh.min_,3)

sections_list = [([1,0,0], dim_min[0], dim_max[0], 10),
                 ([0,1,0], dim_min[1], dim_max[1], 1),
                 ([0,0,1], dim_min[2], dim_max[2], 1)]

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

drawing = dxf.drawing('test.dxf')

for i, section_list in enumerate(section_plane_list):
    for section in section_list:
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
                    drawing.add(dxf.line((y0+x0, z0), (y1+x1, z1), color=7,layer='~plane_{}'.format(i)))
                elif i==1:
                    drawing.add(dxf.line((x0, z0), (x1, z1), color=2,layer='~plane_{}'.format(i)))
                elif i==2:
                    drawing.add(dxf.line((x0, y0), (x1, y1), color=1,layer='~plane_{}'.format(i)))
drawing.save()


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
