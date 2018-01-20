import numpy as np
from stl import mesh
import cncfclib as fc
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from dxfwrite import DXFEngine as dxf



def find3dpoint(p, pf):
    '''function searches for /pf/ vector in an array of vectors /p/
       and returns its top level position index'''

    if len(p.shape)>2:
        pool_A=p[:,0]
        pool_B=p[:,1]
        idx_row0=np.where(np.product(pool_A==pf,axis=1))[0]
        idx_row1=np.where(np.product(pool_B==pf,axis=1))[0]
    else:
        pool_A=p
        idx_row0=np.where(np.product(pool_A==pf,axis=1))[0]
        idx_row1=[]
    return np.hstack([idx_row0, idx_row1]).astype(np.int)


def make_loop(b):
    '''function returns an array of subsequent points from an unsorted
       array of contour sections /b/. An example of the contour sections array:
    p=np.array([[[2, 4,	0], [1,	1, 0]],
                [[1, 2,	0],	[4,	4, 0]],
                [[2, 4, 0], [4, 4, 0]],
                [[1, 1,	0],	[1,	2, 0]]])
    '''
    p=b.astype(np.float)
    p=p.round(4)
    p_out = np.zeros_like(p[:,0])
    idx_row = 0
    idx_col = 0

    p_first = p[idx_row, idx_col]
    p_out[0] = p_first
    for i in np.arange(0, p.shape[0]):
        next_idx_row = np.setdiff1d(find3dpoint(p, p_first), idx_row)[0]
        next_idx_col = find3dpoint(p[next_idx_row], p_first)
        el_next=p[ next_idx_row, next_idx_col^1]
        p_out[i] = el_next
        p_first = el_next
        idx_row = next_idx_row
        idx_col = next_idx_col
    return p_out

def remove_mid_nodes(arr):
    p=arr
    print(p)
    while 1:
        vA=np.roll(p, -1, axis=0)-p
        vB=np.roll(p, -1, axis=0)-np.roll(p, -2, axis=0)
        c=np.linalg.norm(np.cross(vA,vB),axis=1)
        # print(c)
        # print(np.cross(vA,vB))
        idx=np.where(c<2)[0]
        if np.size(idx):
            print(p)
            print(c[idx])
            print(np.size(p[:,0])-1)
            idxi=np.where(idx>=np.size(p[:,0])-1)[0]
            if idxi:
                idx[idxi]=0
            p=np.delete(p,idx+1,0)
        else:
            break
    return p

def make_chains(section_list):
    chain_list
    p_arr =  np.array(section)
    print(p_arr.shape[0])
    print('section: ',i)
    if p_arr.shape[0]>3:
        prof=make_loop(p_arr)
    return chain_list


mesh = mesh.Mesh.from_file('fuselage.stl')
dim_max = mesh.max_
dim_min = mesh.min_

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
    # print('slicing the model')
    #
    for i, (n0, D0) in enumerate(zip(cp_n0_arr, cp_D0_arr)):
        intersect_list = []
        for tri in mesh.vectors:
            #ABC
            P1 = np.vstack(tri).astype(float)
            #CAB
            P2 = np.roll(P1, 1, axis=0)
            # print(ABC)
            # print(BCA)
            intersect = fc.tri_plane_intersect_check(P1, P2, D0, n0)
            if np.size(intersect):
                print(intersect)
                intersect_list.append(intersect)
    # #
        print('profile: {}; sections: {}'.format(i,len(intersect_list)))
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
                    drawing.add(dxf.line((y0+x0, z0), (y1+x1, z1), color=7,layer='plane_{}'.format(i)))
                elif i==1:
                    drawing.add(dxf.line((x0, z0), (x1, z1), color=2,layer='plane_{}'.format(i)))
                elif i==2:
                    drawing.add(dxf.line((x0, y0), (x1, y1), color=1,layer='plane_{}'.format(i)))
drawing.save()


fig = plt.figure()
ax = fig.gca(projection='3d')

for section_list in section_plane_list:
    for i, section in enumerate(section_list):
        p_arr =  np.array(section)
        for row in p_arr:
            if row.shape[0]==2:
                print(row)
                x = row[:,0]
                y = row[:,1]
                z = row[:,2]
                ax.plot(x, y, z,'s-')
plt.show()
