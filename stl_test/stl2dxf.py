import numpy as np
from stl import mesh
import cncfclib as fc
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



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

# B# print(point_plane_dist(P1, D0, n0))
# n_sect=20
# cp_n0_arr = np.array([[1,0,0]]*n_sect)
# cp_D0_arr = np.array([[1,0,0]]*n_sect) *
mesh = mesh.Mesh.from_file('fuselage.stl')
n_sect=2
cp_n0_arr = np.tile([0,0,1],(n_sect,1))
cp_D0_arr = np.tile([0,0,1],(n_sect,1)) * np.linspace(0,250,n_sect)[:,np.newaxis]
# print(cp_D0_arr)
# print(cp_D0_arr)
#
section_list=[]
print('slicing the model')
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

#         if np.size(intersect):
#             # print(intersect)
#             intersect_list.append(intersect)
# #
#     print('profile: {}; sections: {}'.format(i,len(intersect_list)))
#     section_list.append(intersect_list)
#
# # print(section_list[0])
#
# fig = plt.figure()
# ax = fig.gca(projection='3d')
#
# for i, section in enumerate(section_list):
#     p_arr =  np.array(section)
#     # print(p_arr)
#     # print(p_arr.shape[0])
#     # print('section: ',i)
#     if p_arr.shape[0]>3:
#     #     # prof=make_loop(p_arr)
#         # profi=remove_mid_nodes(prof)
#         # print(prof)
#         x = p_arr[:,0,0]
#         y = p_arr[:,0,1]
#         z = p_arr[:,0,2]
#         # print(z)
#         ax.plot(x, y, z,'s-')
#         # x = profi[:,0]
#         # y = profi[:,1]
#         # z = profi[:,2]
#         # ax.plot(x, y, z,'o')
#         # ax.plot(x[[0,-1]], y[[0,-1]], z[[0,-1]], 'o-')
#
#         # test_p_arr=np.vstack(p_arr)
#         # print(test_p_arr)
#         # x = test_p_arr[:,0]
#         # y = test_p_arr[:,1]
#         # z = test_p_arr[:,2]
#         # ax.plot(x, y, z,'x')
# plt.show()
