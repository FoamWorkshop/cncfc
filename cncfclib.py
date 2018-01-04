import sys
import os
import numpy as np
from stl import mesh
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def read_data(f_name, msg='False'):

    data = []

    if os.path.isfile(f_name):

        with open(f_name, 'r') as f:
            for line in f:
                tmp = line.split()
                x, y, z = 0, 0, 0

                if len(tmp) == 1:
                    x = tmp[0]
                    y = 0
                    z = 0
                elif len(tmp) == 2:
                    x, y = tmp
                    z = 0
                elif len(tmp) == 3:
                    x, y, z = tmp

                data.append([float(x), float(y), float(z)])
        if msg:
            print("{0:<24} -> {1} knots".format(f_name, len(data)))

    else:
        if msg:
            print('{0} not found'.format(f_name))

    return data


def write_data(f_name, data, msg='False'):
    i = 0
    print('assigned file name', f_name)
    # if os.path.isfile(f_name):
    #
    #     bak_f_name = f_name + '.00.bak'
    #
    #     while os.path.isfile(bak_f_name):
    #         bak_f_name = '{0}.{1:{fill}>2s}.bak'.format(
    #             f_name, str(i), fill='0')
    #         i += 1
    #
    #     os.rename(f_name, bak_f_name)
    #
    #     if msg:
    #         print("{0:<24} -> {1}".format(f_name, bak_f_name))

    with open(f_name, 'w') as f:
        for line in data:
            x, y, z = line
            f.write('{0:.6f} {1:.6f} {2:.6f}\n'.format(
                float(x), float(y), float(z)))

    if msg:
        print("{0:<24} <- {1} knots".format(f_name, len(data)))



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

# def make_loop(b):
#     '''function returns an array of subsequent points from an unsorted
#        array of contour sections /b/. An example of the contour sections array:
#     p=np.array([[[2, 4,	0], [1,	1, 0]],
#                 [[1, 2,	0],	[4,	4, 0]],
#                 [[2, 4, 0], [4, 4, 0]],
#                 [[1, 1,	0],	[1,	2, 0]]])
#     '''
#     p=b.astype(np.float)
#     p=p.round(4)
#     p_out = np.zeros_like(p[:,0])
#     idx_row = 0
#     idx_col = 0
#
#     p_first = p[idx_row, idx_col]
#     p_out[0] = p_first
#     for i in np.arange(0, p.shape[0]):
#         next_idx_row = np.setdiff1d(find3dpoint(p, p_first), idx_row)[0]
#         next_idx_col = find3dpoint(p[next_idx_row], p_first)
#         el_next=p[ next_idx_row, next_idx_col^1]
#         p_out[i] = el_next
#         p_first = el_next
#         idx_row = next_idx_row
#         idx_col = next_idx_col
#     return p_out

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
    if np.size(p):
        cross_vect=np.apply_along_axis(line_plane_intersect,1,L[p,:],D0,n0)
    return cross_vect

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
def cartesian2cylyndrical(data_points,sections):
    # p=p.round(2)
    p=np.unique(data_points,axis=0)
    p_size=np.size(p)

    x = p[:,0]
    y = p[:,1]
    z = p[:,2]

    r =  np.linalg.norm(np.vstack((x,y)),axis=0)
    th = np.arctan2(x,y)
    pos = np.vstack([r, th, z]).T
    ind = np.lexsort((th,z))
    pos=pos[ind].reshape((-1,sections-1,3))
    # print(np.shape(pos))
    for i in np.arange(np.shape(pos)[0]):
        # print(pos[i,:,1])
        ind=np.argsort(pos[i,:,1])
        pos[i]=pos[i,ind,:]
        # print('section ',i)
        # print(pos[i])
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

    return np.vstack([rho * np.cos(phi),
                      rho * np.sin(phi),
                      z]).T

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

def plot_loft_paths(data):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for i in np.arange(np.shape(data)[1]):
        # ax.plot(spars[:,0],spars[:,1],spars[:,2],'x-')
        ax.plot(data[:,i,0],data[:,i,1],data[:,i,2],'x-')
    plt.show()
