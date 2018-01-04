import numpy as np
from stl import mesh
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import cncfclib as cf
from scipy.interpolate import griddata

def v2v_dist(v1, v2):
    return np.linalg.norm(v2-v1)

def point_plane_dist(P1, D0, n0):
    v1 = P1-D0
    n_norm = n0/np.linalg.norm(n0)
    return np.dot(n_norm,v1)
    # +D0

def line_plane_intersect_d(L0, L1, D0, n0):
    l_norm = (L1-L0)
    n_norm = n0/np.linalg.norm(n0)
    d = np.dot((D0 - L0),n_norm)/np.dot(l_norm, n_norm)
    return d

def line_plane_intersect(L, D0, n0):
    L0 = L[:3]
    L1 = L[3:]
    l = (L1-L0)
    n_norm = n0/np.linalg.norm(n0)
    d = np.dot((D0 - L0),n_norm)/np.dot(l, n_norm)
    return L0 + d * l
def cartesian2cylyndrical(data_points):
    # p=p.round(2)
    p=np.unique(data_points,axis=0)
    p_size=np.size(p)
    r = np.linalg.norm(p[:,:2],axis=1)
    th = np.arctan2(p[:,0],p[:,1])
    z = p[:,2]
    pos = np.vstack([r, th, z]).T
    ind = np.lexsort((pos[:,1],pos[:,2]))
    pos=pos[ind].reshape((-1,16,3))
    print(np.shape(pos))
    for i in np.arange(np.shape(pos)[0]):
        # print(pos[i,:,1])
        ind=np.argsort(pos[i,:,1])
        pos[i]=pos[i,ind,:]
        print('section ',i)
        print(pos[i])
        # print(pos)
    return pos

def cylyndrical2cartesian(p):
    if p.size>3:
        rho=p[:,0]
        phi=p[:,1]
        z=p[:,2]
    else:
        rho=p[0]
        phi=p[1]
        z=p[2]

    return np.vstack([rho * np.cos(phi),rho * np.sin(phi),z]).T

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

def find_nearest(p, pf):
    idx=1
    return idx

def round_arr(p, dec):
    return np.apply_along_axis(np.round,0,p,dec)

def slice_mesh(mesh, cp_n0_arr, cp_D0_arr):
    section_list=[]
    print('slicing the model')
    for i, (n0, D0) in enumerate(zip(cp_n0_arr, cp_D0_arr)):
        intersect_list = []
        for tri in mesh.vectors:
            ABC = np.vstack(tri)
            L=np.hstack([ABC,np.roll(ABC, -1, axis=0)]).astype(float)
            intersect=tri_plane_intersect_check(L,D0,n0)
            if np.size(intersect):
                intersect_list.append(intersect)

        print('profile: {}; sections: {}'.format(i,len(intersect_list)))
        section_list.append(intersect_list)
    return section_list

def make_loop(b):
    '''function returns an array of subsequent points from an unsorted
       array of contour sections /b/. An example of the contour sections array:
    p=np.array([[[2, 4,	0], [1,	1, 0]],
                [[1, 2,	0],	[4,	4, 0]],
                [[2, 4, 0], [4, 4, 0]],
                [[1, 1,	0],	[1,	2, 0]]])
    '''
    p_out=[]
    if b.size:
        print(b.size)
        # print(b)
        p=b.astype(np.float)
        p=b.round(4)
        p_out = np.zeros_like(p[:,0])
        idx_row = 0
        idx_col = 0

        p_first = p[idx_row, idx_col]
        p_out[0] = p_first
        for i in np.arange(0, p.shape[0]):
            next_idx_row = np.setdiff1d(find3dpoint(p, p_first), idx_row)[0]
            next_idx_col = find3dpoint(p[next_idx_row], p_first)
            el_next=p[ next_idx_row, next_idx_col^1]
            print(el_next)
            p_out[i] = el_next
            p_first = el_next
            idx_row = next_idx_row
            idx_col = next_idx_col
    return p_out

def plot_section(section):
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    for edge in section:
        p_arr =  np.vstack(edge)

        x = p_arr[:,0]
        y = p_arr[:,1]
        z = p_arr[:,2]
        ax.plot(x, y, z)
            # ax.plot(x[[0,-1]], y[[0,-1]], z[[0,-1]], 'o-')
    plt.show()

def tri_plane_intersect_check(L,D0,n0):
    P1 = L[:,:3]
    P2 = L[:,3:]
    P1_dist=np.apply_along_axis(point_plane_dist,1,P1,D0,n0)
    P2_dist=np.apply_along_axis(point_plane_dist,1,P2,D0,n0)
    p=np.where(P1_dist*P2_dist<=0)[0]
    cross_vect=np.array([])
    if p.size:
        cross_vect=np.apply_along_axis(line_plane_intersect,1,L[p,:],D0,n0)
    return cross_vect

def make_chains(section_list):
    chain_list
    p_arr =  np.array(section)
    print(p_arr.shape[0])
    print('section: ',i)
    if p_arr.shape[0]>3:
        prof=make_loop(p_arr)
    return chain_list


mesh = mesh.Mesh.from_file('fuselage_rot.stl')

# sect_n=10
# sect_space = np.linspace(0.001,200,sect_n)
# cp_n0_Barr = np.array([[0,0,1]]*sect_n)
# cp_D0_arr = np.array([[0,0,1]]*sect_n) * np.vstack(sect_space)
# section_list=slice_mesh(mesh, cp_n0_arr, cp_D0_arr)

sect_n=1
sect_space = np.linspace(0,3,sect_n)
cp_n0_arr = np.array([[0,1,0]]*sect_n)
cp_D0_arr = np.array([[0,1,0]]*sect_n) * np.vstack(sect_space)
# section_sym = slice_mesh(mesh, cp_n0_arr, cp_D0_arr)
v_arr = np.round(np.vstack(mesh.vectors).astype(float), decimals=1)
pos = cartesian2cylyndrical(v_arr)
# print(v_arr_cyl)
# th = v_arr_cyl[:,1]
# z = v_arr_cyl[:,2]
# r = v_arr_cyl[:,0]


# print(v_arr_cyl)
# th=np.linspace(-np.pi*0.75,np.pi*0.75,36)
# z = 30*np.ones_like(th)
# x = np.vstack([th,z]).T
# r = griddata( v_arr_cyl[:,1:] , v_arr_cyl[:,0], x, method='linear')
# print(r)
# print(np.argmax(v_arr_cyl[:,1]))
# print(v_arr_cyl[np.argmax(v_arr_cyl[:,1])])
# print(np.pi)

# fig = plt.figure()
# ax = fig.gca()
# for i in np.arange(np.shape(pos)[0]):
#     ax.plot(pos[i,:,1],pos[i,:,2],'x')
# plt.show()
#
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# for i in np.arange(np.shape(pos)[0]):
#     spars = cylyndrical2cartesian(pos[i])
#     ax.plot(spars[:,0],spars[:,1],spars[:,2],'x-')
# plt.show()

fig = plt.figure()
ax = fig.gca(projection='3d')
for i in np.arange(np.shape(pos)[1]):
    spars = cylyndrical2cartesian(pos[:,i])
    ax.plot(spars[:,0],spars[:,1],spars[:,2],'x-')
plt.show()


# fig = plt.figure()
# ax = fig.gca()
# # ax = fig.gca(projection='3d')
# for z_mem in np.linspace(0,300,20):
#     the = np.linspace(-np.pi*0.6, np.pi*0.6,17)
#     z = np.ones_like(the)*z_mem
#     x=np.vstack((the,z)).T
#     r = griddata( v_arr_cyl[:,1:] , v_arr_cyl[:,0], x, method='nearest')
#     v_arr_car = cylyndrical2cartesian(np.vstack([r,the,z]).T)
#     # print(r)
#     # ax.plot(r, the,'o-')
#     ax.plot(v_arr_cyl[:,2],v_arr_cyl[:,1],'x')
#     # plt.ylim(-np.pi,np.pi)
#     # plt.xlim(0,15)
#     # plt.grid(True)
#     # ax.plot(v_arr_cyl[:,0], v_arr_cyl[:,1], v_arr_cyl[:,2],'o')
#     # ax.plot(r, the, z,'s-')
#     # ax.plot(v_arr_car[:,0], v_arr_car[:,1], v_arr_car[:,2],'o-')
# plt.show()

    # v_arr_car = cylyndrical2cartesian(grid)
    # print(v_arr_car)
# print('r: ',r)
# # print(r.T)
# v_arr_car = cylyndrical2cartesian(np.hstack([np.vstack(r),x]))
#
# print(v_arr_car)
#
# # print(v_arr_cyl)
# # print(v_arr.shape)
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# # v_arr_cyl
# x = v_arr_cyl[:,0]
# y = v_arr_cyl[:,1]
# z = v_arr_cyl[:,2]
# # ax.plot(x, y, z)
# # x = v_arr[:,0]
# # y = v_arr[:,1]
# # z = v_arr[:,2]
# ax.plot(x, y, z,'o-')
# plt.show()



# plot_section(section_list[0])

# fig = plt.figure()
# ax = fig.gca(projection='3d')

# for i, section in enumerate(section_list):
#     p_arr =  np.array(section)
#     print(p_arr.shape[0])
#     print('section: ',i)
#     if p_arr.shape[0]>3:
#         prof=make_loop(p_arr)
#         x = prof[:,0]
#         y = prof[:,1]
#         z = prof[:,2]
#         ax.plot(x, y, z)
#         ax.plot(x[[0,-1]], y[[0,-1]], z[[0,-1]], 'o-')
# for i, section in enumerate(section_sym):
#     p_arr =  np.array(section)
#     print(p_arr.shape[0])
#     print('section: ',i)
#     if p_arr.shape[0]>3:
#         prof=make_loop(p_arr)
#         x = prof[:,0]
#         y = prof[:,1]
#         z = prof[:,2]
#         ax.plot(x, y, z)
#         ax.plot(x[[0,-1]], y[[0,-1]], z[[0,-1]], 'o-')
#
#
# plt.show()
