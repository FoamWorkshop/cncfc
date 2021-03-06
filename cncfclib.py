#!/usr/bin/python3
import sys
import os
import numpy as np
from stl import mesh
import dxfgrabber
import ezdxf
import pprint as pp
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from collections import Iterable
from scipy.linalg import expm, norm
from scipy import spatial
from itertools import groupby
from operator import itemgetter
import collections
from collections import defaultdict
from collections import OrderedDict
import re
from copy import deepcopy
# def flatten(items):
#     """Yield items from any nested iterable; see REF."""
#     for x in items:
#         if isinstance(x, Iterable) and not isinstance(x, (str, bytes)):
#             yield from flatten(x)
#         else:
#             yield x

def nested_dict():
    return defaultdict(nested_dict)

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

    with open(f_name, 'w') as f:
        for line in data:
            x, y, z = line
            f.write('{0:.6f} {1:.6f} {2:.6f}\n'.format(
                float(x), float(y), float(z)))

    if msg:
        print("{0:<24} <- {1} knots".format(f_name, len(data)))

def knots2gcode(ct_pathxy, ct_pathuv, ct_path_r, name='gcode', global_header='False', subset_header='False'):
    f_name = name + '.ngc'
    with open(f_name, 'w') as f:
        if subset_header:
            f.write("o<{0}> sub\n".format(name))

        if global_header:

            f.write("G21 (Units in millimeters)\n")
            f.write("G90 (Absolute programming)\n")
            f.write("G40 (Cancel radius comp.)\n")
            f.write("G49 (Cancel length comp.)\n")
            f.write("F200 (Feed rate)\n")
    ##############################################
        f.write("\n(-----CT PATH-from: {0}-----)\n".format(name))
        if ct_path_r:
            f.write('G1 B{0:8.3f}\n'.format(ct_path_r[0][0]))
            for varxy, varuv, var_r in zip(ct_pathxy, ct_pathuv, ct_path_r):
                f.write('G1 X{0:8.3f} Y{1:8.3f} U{2:8.3f} V{3:8.3f} B{4:8.3f}\n'.format(
                    varxy[0], varxy[1], varuv[0], varuv[1], var_r[0]))
        else:
            for varxy, varuv in zip(ct_pathxy, ct_pathuv):
                f.write('G1 X{0:8.3f} Y{1:8.3f} U{2:8.3f} V{3:8.3f}\n'.format(
                    varxy[0], varxy[1], varuv[0], varuv[1]))
    ##############################################
        if subset_header:
            f.write("\no<{0}> endsub\n".format(name))
            f.write("M2 (program end)")

        if global_header:
            f.write("M2 (program end)")

def savegcode(data, name='gcode', global_header=False, subset_header=False):

    f_name = name + '.ngc'
    with open(f_name, 'w') as f:
        if subset_header:
            f.write("o<{0}> sub\n".format(name))

        if global_header:

            f.write("G21 (Units in millimeters)\n")
            f.write("G90 (Absolute programming)\n")
            f.write("G40 (Cancel radius comp.)\n")
            f.write("G49 (Cancel length comp.)\n")
            f.write("F200 (Feed rate)\n")
    ##############################################
        f.writelines(flatten(data))
        # f.writelines(data[3])
    ##############################################
        if subset_header:
            f.write("\no<{0}> endsub\n".format(name))
            f.write("M2 (program end)")

        if global_header:
            f.write("M2 (program end)")


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
    if np.linalg.norm(n0)==0:
        print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
    n_norm = n0 / np.linalg.norm(n0)
    d = np.dot((D0 - L0),n_norm) / np.dot(l, n_norm)
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


def tri_plane_intersect_check(P1,P2,D0,n0):
    # P1 = L[:,:3]
    # P2 = L[:,3:]
    P1_dist=np.apply_along_axis(point_plane_dist,1,P1,D0,n0)
    P2_dist=np.apply_along_axis(point_plane_dist,1,P2,D0,n0)
    p_idx=np.where(P1_dist*P2_dist<=0)[0]
    #
    if np.size(p_idx) :
        # for idx in p_idx:
            # line_plane_intersect(P1, P2, D0, dn)
        cross_vect=np.apply_along_axis(line_plane_intersect,1,np.column_stack((P1,P2))[p_idx],D0,n0)
        ###rid of nans!!!!!
        valid_idx = np.where(np.isfinite(cross_vect[:,0]))[0]
        # if np.size(valid_idx):
            # print('00000000000000000000000000000',valid_idx)
            # print(cross_vect)
        cv = cross_vect[valid_idx]
    else:
        cv=np.array([])
    return cv

def line_plane_intersect(L0, L1, D0, n0):
    """Calculate 3D line-plane intersection point
    Parameters:
        L0, L1 - points on the line
        D0 - point on the plane
        n0 - plane normal
    Returns:
        P - intersetion point. If the line is parallel to the plane,
        function returns [inf, inf, inf]
    """
    l = (L1 - L0)
    n_norm = n0/np.linalg.norm(n0)
    dot_prod = np.dot(l, n_norm)

    if dot_prod == 0:
        P = np.array([np.inf,np.inf,np.inf])
    else:
        d = np.dot((D0 - L0),n_norm)/dot_prod
        P = L0 + d * l
    return P

def tri_plane_intersect(tri0, D0, n0):
    """Calculate 3D triangle-plane intersection segment
    Parameters:
        tri0 - array of the triangle vartice coordinates
            [[ax,ay,az],[bx,by,bz],[cx,cy,cz]]
        D0 - point on the plane
        n0 - plane normal

    Returns:
        [S1, S2] - a list of intersection point coordinates (arrays)
        defining the intersection segment, if no intersection occures,
        an empty list is returned
    """
    tri1 = np.roll(tri0, 1, axis=0)
    tri0_dist=np.apply_along_axis(point_plane_dist,1,tri0,D0,n0)
    tri1_dist=np.apply_along_axis(point_plane_dist,1,tri1,D0,n0)
    idx = np.where(tri0_dist*tri1_dist<=0)[0]
    #
    S_list=[]
    if np.size(idx):
        for L0, L1 in zip(tri0[idx], tri1[idx]):
            S = line_plane_intersect(L0, L1, D0, n0)

            if np.isinf(S[0]):
                S_list = [L0,L1]
                break
            else:
                S_list.append(S)

    return S_list



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

# def line_plane_intersect(L, D0, n0):
#     L0 = L[:3]
#     L1 = L[3:]
#     l = (L1 - L0)
#     n_norm = n0/np.linalg.norm(n0)
#     dot_prod = np.dot(l, n_norm)
#
#     if dot_prod !=0:
#         d = np.dot((D0 - L0),n_norm)/dot_prod
#         res = L0 + d * l
#     else:
#         res = np.full((1,3),np.inf)
#     return res



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
            # print(el_next)
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

# def tri_plane_intersect_check(L,D0,n0):
#     P1 = L[:,:3]
#     P2 = L[:,3:]
#     P1_dist=np.apply_along_axis(point_plane_dist,1,P1,D0,n0)
#     P2_dist=np.apply_along_axis(point_plane_dist,1,P2,D0,n0)
#     p=np.where(P1_dist*P2_dist<=0)[0]
#     cross_vect=np.array([])
#     if p.size:
#         cross_vect=np.apply_along_axis(line_plane_intersect,1,L[p,:],D0,n0)
#     return cross_vect

def plot_loft_paths(data):
    fig = plt.figure()
    # plt.ion()
    ax = fig.gca(projection='3d')
    for i in np.arange(np.shape(data)[1]):
        # ax.plot(spars[:,0],spars[:,1],spars[:,2],'x-')
        ax.plot(data[:,i,0],data[:,i,1],data[:,i,2],'x-')
    plt.show()

def plot_surf(x, y, z):
    fig = plt.figure()
    # plt.ion()
    ax = fig.gca(projection='3d')
    for s1, s2, s3 in zip(x,y,z):
        # ax.plot(spars[:,0],spars[:,1],spars[:,2],'x-')
        ax.plot(s1, s2, s3,'x-')
    plt.show()

def angl_conv(x):
    out = x

    if x<0:
        out = 2*pi + x

    return out

# vangl=np.vectorize(angl_conv)
def dist(P1, P2):
    return np.sqrt(np.sum(((P1-P2)**2)))

def radius_segment(P1,P2,P3):
    P4 = np.abs( (P2[:,1]-P1[:,1]) * P3[:,0] - (P2[:,0]-P1[:,0]) * P3[:,1] + P2[:,0]*P1[:,1] - P2[:,1]*P1[:,0]) / np.sqrt( (P2[:,1] - P1[:,1])**2 + (P2[:,0]-P1[:,0])**2 )
    return np.vstack(P4)

def radius(P1,P2):
    r= np.sqrt((P1[:,0] - P2[:,0])**2 + (P1[:,1] - P2[:,1])**2)
    return np.vstack(r)

def angle_segment(P1,P2,P3):
    '''   P2
         /
        /
       P4.
      /       .
    P1             P3
    '''
    k = ((P2[:,1]-P1[:,1]) * (P3[:,0]-P1[:,0]) - (P2[:,0]-P1[:,0]) * (P3[:,1]-P1[:,1])) / ((P2[:,1]-P1[:,1])**2 + (P2[:,0]-P1[:,0])**2)
    P4 = np.vstack([P3[:,0] - k * (P2[:,1]-P1[:,1]), P3[:,1] + k * (P2[:,0]-P1[:,0])]).T
    angl = np.arctan2(P4[:,1] - P3[:,1], P4[:,0] - P3[:,0])
    angl = np.vstack(angl)
    return angl#np.apply_along_axis(lambda x:x if x>=0 else 2*pi+x, 1, angl)

def cross_point(P1,P2,P3):
    '''   P2
         /
        /
       P4.
      /       .
    P1             P3
    '''
    k = ((P2[:,1]-P1[:,1]) * (P3[:,0]-P1[:,0]) - (P2[:,0]-P1[:,0]) * (P3[:,1]-P1[:,1])) / ((P2[:,1]-P1[:,1])**2 + (P2[:,0]-P1[:,0])**2)
    P4 = np.vstack([P3[:,0] - k * (P2[:,1]-P1[:,1]), P3[:,1] + k * (P2[:,0]-P1[:,0])]).T


    return P4#np.apply_along_axis(lambda x:x if x>=0 else 2*pi+x, 1, angl)

def angle_test(P1,P2,P3):
    '''   P3
         /
        / angl
       P2-------P1

    v1 = P1-P2 <- ref. vector
    v2 = P3-P2
    '''
    v1 = P1-P2
    v2 = P3-P2

    angl1 = np.vstack(np.arctan2(v1[:,1],v1[:,0]))
    angl2 = np.vstack(np.arctan2(v2[:,1],v2[:,0]))
    dangl = np.apply_along_axis(lambda x:x if x>=0 else 2*pi+x, 1, angl2-angl1)

    return np.vstack(dangl)

def angle_atan2(P1,P2,P3):
    '''   P3
         /
        / angl
       P2-------P1

    v1 = P1-P2 <- ref. vector
    v2 = P3-P2
    '''
    v1 = P1-P2
    v2 = P3-P2

    angl1 = np.vstack(np.arctan2(v1[:,1],v1[:,0]))
    angl2 = np.vstack(np.arctan2(v2[:,1],v2[:,0]))
    dangl = np.apply_along_axis(lambda x:x if x>=0 else 2*pi+x, 1, angl2-angl1)
    dangl = np.apply_along_axis(lambda x:x if x<=pi else x-2*pi, 1, dangl)

    return np.vstack(dangl)

def coords2file(name, coords_XU, coords_YV):
    pref_1='xyuv_'
    pref_2='r_'

    f = open(name, 'w')
    for XU, YV in zip(coords_XU,coords_YV):
        f.write('{0:.3f} {1:.3f}\n'.format(XU, YV))

    f.close()

def Rcoords2file(name,coords_R):
    pref_1='xyuv_'
    pref_2='r_'

    f = open(name, 'w')
    for R in coords_R:
        f.write('{0:.3f}\n'.format(R))
    f.close()

def coords2gcode(x_arr,y_arr,u_arr,v_arr,b_arr):
    # pref_1='xyuv_'
    # pref_2='r_'
    # if !b:
    #     b_arr = np.zeros_like(x)
    gcode=[]
    for row in np.column_stack((x_arr, y_arr, u_arr, v_arr, b_arr)):
        gcode.append('G1 X{0[0]:8.3f} Y{0[1]:8.3f} U{0[2]:8.3f} V{0[3]:8.3f} B{0[4]:8.3f}\n'.format(row))
    return gcode


def list_entities(dxf):
    dxf_summary = [shape.dxftype for shape in dxf.entities]
    print('{0:<10}: {1}'.format('LINES', dxf_summary.count('LINE')))
    print('{0:<10}: {1}'.format('ARCS', dxf_summary.count('ARC')))

def knots_dict(knots_list):
    return [[i, var] for i, var in enumerate(list(set(knots_list)))]

def interp_points(seg_P1_X_list, seg_P1_Y_list, seg_O_X_list, seg_O_Y_list, C_Z, n_sect, interp_meth, name, layer_list):
    P1_X=np.vstack(seg_P1_X_list)
    P1_Y=np.vstack(seg_P1_Y_list)
    O_X=np.vstack(seg_O_X_list)
    O_Y=np.vstack(seg_O_Y_list)

    yv = np.hstack((np.linspace(C_Z[0],C_Z[-1],n_sect),C_Z))
    yv = np.unique(yv)
    yv = np.sort(yv)
    #
    n_spokes = len(seg_P1_X_list[0])
    n_sections = len(yv)

    truss_list=[]    # print P1_X

    interp_P1_X_list = []
    interp_P1_Y_list = []
    interp_O_X_list = []
    interp_O_Y_list = []
    #
    for i in range(n_spokes):

        interp_P1_X = interpolate.interp1d(C_Z, P1_X[:,i],kind=interp_meth)(yv)
        interp_P1_Y = interpolate.interp1d(C_Z, P1_Y[:,i],kind=interp_meth)(yv)
        interp_O_X =  interpolate.interp1d(C_Z, O_X[:,0],kind=interp_meth)(yv)
        interp_O_Y =  interpolate.interp1d(C_Z, O_Y[:,0],kind=interp_meth)(yv)
    #
        interp_P1_X_list.append(interp_P1_X)
        interp_P1_Y_list.append(interp_P1_Y)
        interp_O_X_list.append(interp_O_X)
        interp_O_Y_list.append(interp_O_Y)

    #print np.array([interp_O_X,interp_O_Y]).T
    #
    interp_P2_X_list= interp_P1_X_list[1:] + interp_P1_X_list[:1]
    interp_P2_Y_list= interp_P1_Y_list[1:] + interp_P1_Y_list[:1]

    P1_X = np.vstack([interp_P1_X_list]).T
    P1_Y = np.vstack([interp_P1_Y_list]).T
    P2_X = np.vstack([interp_P2_X_list]).T
    P2_Y = np.vstack([interp_P2_Y_list]).T
    O_X = np.vstack([interp_O_X_list]).T
    O_Y = np.vstack([interp_O_Y_list]).T

    C_a = np.ones([n_sections,n_spokes])
    C_r = np.ones([n_sections,n_spokes])

    P1 = np.vstack([P1_X[:,0], P1_Y[:,0]]).T
    P2 = np.vstack([P2_X[:,0], P2_Y[:,0]]).T
    O  = np.vstack([O_X[:,0], O_Y[:,0]]).T

    P_ref = np.ones((len(yv) , 2)) * np.array([[1,0]])
    C_S = cross_point(P1, P2, O)
    C_a_ref = angle_test( P_ref, O, C_S) #* 180/pi

    for i in range(n_spokes):
        # print i
        P1 = np.vstack([P1_X[:,i], P1_Y[:,i]]).T
        P2 = np.vstack([P2_X[:,i], P2_Y[:,i]]).T
        O = np.vstack([O_X[:,i], O_Y[:,i]]).T
        C = cross_point(P1, P2, O)
        C_r[:,i]=radius(C,O)[:,0]
        C_a[:,i]=(angle_test( C_S, O, C) + C_a_ref)[:,0]* 180/pi

    # print C_a
        # C_a[:,i] = C_a[:,i]#+C_a_ref
    # print C_r
    # print C_a
    # # #
    cut_in_swing = True
    for i in range(n_spokes):

        f_name_C_r1='{0}_xyuv_{1:{fill}{align}4}.knt'.format(layer_list,i,fill='0',align='>')
        f_name_C_a1='{0}_b_{1:{fill}{align}4}.knt'.format(layer_list,i,fill='0',align='>')

        if cut_in_swing and i%2:
            coords2file(f_name_C_r1, np.flipud(C_r[:,i]), np.flipud(yv))
            Rcoords2file(f_name_C_a1, np.flipud(C_a[:,i]))
        else:
            coords2file(f_name_C_r1, C_r[:,i], yv)
            Rcoords2file(f_name_C_a1, C_a[:,i])

def interp_points_cyl(strokes, name, layer_list):
# def interp_points(seg_P1_X_list, seg_P1_Y_list, seg_O_X_list, seg_O_Y_list, C_Z, n_sect, interp_meth, name, layer_list):
    # P1_X=np.vstack(seg_P1_X_list)
    # P1_Y=np.vstack(seg_P1_Y_list)
    # O_X=np.vstack(seg_O_X_list)
    # O_Y=np.vstack(seg_O_Y_list)
    #
    # yv = np.hstack((np.linspace(C_Z[0],C_Z[-1],n_sect),C_Z))
    # yv = np.unique(yv)
    # yv = np.sort(yv)
    # #
    # n_spokes = len(seg_P1_X_list[0])
    # n_sections = len(yv)
    #
    # truss_list=[]    # print P1_X
    #
    # interp_P1_X_list = []
    # interp_P1_Y_list = []
    # interp_O_X_list = []
    # interp_O_Y_list = []
    # #
    # for i in range(n_spokes):
    #
    #     interp_P1_X = interpolate.interp1d(C_Z, P1_X[:,i],kind=interp_meth)(yv)
    #     interp_P1_Y = interpolate.interp1d(C_Z, P1_Y[:,i],kind=interp_meth)(yv)
    #     interp_O_X =  interpolate.interp1d(C_Z, O_X[:,0],kind=interp_meth)(yv)
    #     interp_O_Y =  interpolate.interp1d(C_Z, O_Y[:,0],kind=interp_meth)(yv)
    # #
    #     interp_P1_X_list.append(interp_P1_X)
    #     interp_P1_Y_list.append(interp_P1_Y)
    #     interp_O_X_list.append(interp_O_X)
    #     interp_O_Y_list.append(interp_O_Y)
    #
    # #print np.array([interp_O_X,interp_O_Y]).T
    # #
    # interp_P2_X_list= interp_P1_X_list[1:] + interp_P1_X_list[:1]
    # interp_P2_Y_list= interp_P1_Y_list[1:] + interp_P1_Y_list[:1]
    #
    # P1_X = np.vstack([interp_P1_X_list]).T
    # P1_Y = np.vstack([interp_P1_Y_list]).T
    # P2_X = np.vstack([interp_P2_X_list]).T
    # P2_Y = np.vstack([interp_P2_Y_list]).T
    # O_X = np.vstack([interp_O_X_list]).T
    # O_Y = np.vstack([interp_O_Y_list]).T
    #
    # C_a = np.ones([n_sections,n_spokes])
    # C_r = np.ones([n_sections,n_spokes])

    # for i in np.arange(np.shape(data)[1]):
    #     # ax.plot(spars[:,0],spars[:,1],spars[:,2],'x-')
    #     ax.plot(data[:,i,0],data[:,i,1],data[:,i,2],'x-')
    # print(strokes)
    S1 = strokes[:,:,:2]
    S2 = np.roll(strokes,-1, axis=0)
    # O = np.zeros_like(P1)
    #
    # print(P1, P2, O)
    #
    # P_ref = P1[0]
    # C_S = cross_point(P1[0], P2[0], O[0])
    # print(C_S)
    # print(P_ref)
    # print(O)
    # C_a_ref = angle_test( P_ref, O[0], C_S)
    #
    #
    # C_a = np.zeros(np.shape(strokes)[:2])
    # C_r = np.zeros(np.shape(strokes)[:2]).T
    #
    # for i in range(np.shape(P1)[0]):
    #     C = cross_point(P1[i], P2[i], O[i])
    #     # print(C)
    #     print(radius(C, O[i]))
    #     C_r[:,i] = radius(C, O[i]).T
    #     # print('CADASDASDA',C_a)
    #     # C_a[i]=(angle_test( C_S, O[i], C[i]) + C_a_ref)[:,0] * 180/np.pi
    # print(C_r)
    #
    # # print C_a
    #     # C_a[:,i] = C_a[:,i]#+C_a_ref
    # # print C_a
    # # print C_r
    # # #
    # cut_in_swing = True
    # for i in range(n_spokes):
    #
    #     f_name_C_r1='{0}_xyuv_{1:{fill}{align}4}.knt'.format(layer_list,i,fill='0',align='>')
    #     f_name_C_a1='{0}_b_{1:{fill}{align}4}.knt'.format(layer_list,i,fill='0',align='>')
    #
    #     if cut_in_swing and i%2:
    #         coords2file(f_name_C_r1, np.flipud(C_r[:,i]), np.flipud(yv))
    #         Rcoords2file(f_name_C_a1, np.flipud(C_a[:,i]))
    #     else:
    #         coords2file(f_name_C_r1, C_r[:,i], yv)
    #         Rcoords2file(f_name_C_a1, C_a[:,i])

def interp_points_fc(seg_P1_X_list, seg_P1_Y_list, C_Z, n_sect, interp_meth, name):
    print('generate freecad file')
    n_spokes = len(seg_P1_X_list[0])
    print(n_spokes)

    P1_X=np.vstack(seg_P1_X_list)
    P1_Y=np.vstack(seg_P1_Y_list)
    yv = np.hstack((np.linspace(C_Z[0],C_Z[-1],n_sect),C_Z))
    yv = np.unique(yv)
    yv = np.sort(yv)


    truss_list=[]
    # print C_Z
    for i in range(n_spokes):


        if interp_meth=='spline':
            path_P1_X = CubicSpline(C_Z, P1_X[:,i])(yv)
            path_P1_Y = CubicSpline(C_Z, P1_Y[:,i])(yv)

        else:
            path_P1_X = interpolate.interp1d(C_Z, P1_X[:,i],kind=interp_meth)(yv)
            path_P1_Y = interpolate.interp1d(C_Z, P1_Y[:,i],kind=interp_meth)(yv)


        truss_list.append(np.vstack([path_P1_X,path_P1_Y,yv]))

    doc = FreeCAD.newDocument()
    myPart = doc.addObject('Part::Feature', 'trus')

    wire_list_1=[]
    for truss in truss_list:
        points_1 = np.transpose(truss)
        line_list = []

        for p1, p2 in zip(list(points_1[:-1,:]),list(points_1[1:,:])):
            # print(tuple(p1),tuple(p2))
            line_list.append(Part.makeLine(tuple(p1), tuple(p2)))


        wire_list_1.append(Part.Wire(line_list))
    wire_list_2= wire_list_1[-1:] + wire_list_1[:-1]

    for w1, w2 in zip(wire_list_1, wire_list_2):
        myPart = doc.addObject('Part::Feature', 'truss')
        myPart.Shape = Part.makeRuledSurface(w1,w2)
    doc.saveAs(name +'.fcstd')

def angle(v, v_ref, norm = True):
    ''' normalized angle (0..2pi) between 2 vectors\n
    input: v - ref. vector, v_ref - second vector:\n
         v2
         /
        / angl
       -------> v1
    '''
    angl = np.arctan2(v_ref[1], v_ref[0]) - np.arctan2(v[1], v[0])

    if norm:
        if angl < 0:
            angl += 2 * np.pi

    return angl

def cross_point(P1, P2, Q1, Q2):
    '''attachement point, vector and distance of a segment /z/ perpenticular to 2 vectors\n
    input: P1, P2 - u vector points, Q1, Q2 - v vector points\n
    function returns a tupple: PS attachment point on p line, PSQS vector, length of the segment /z/
    ref. geomalgorithms.com #distance between 3d lines & segments#
    '''
    u=P2-P1
    v=Q2-Q1
    u=u/np.linalg.norm(u)
    v=v/np.linalg.norm(v)
    w0=P1-Q1
    a=np.dot(u,u)
    b=np.dot(u,v)
    c=np.dot(v,v)
    d=np.dot(u,w0)
    e=np.dot(v,w0)
    sc=((b*e)-(c*d))/((a*c)-b**2)
    tc=((a*e)-(b*d))/((a*c)-b**2)

    PS = P1 + sc * u
    QS = Q1 + tc * v
    slope_v = QS - PS
    d = np.linalg.norm(slope_v)
    return (PS, slope_v, d)

def normv(v):
    return v/np.linalg.norm(v)

def proj_vector2plane(u, n, y_dir = np.array([0,0,1])):
    '''vector to plane projection.
    input: u - vector, n - plane normal, O - attachement point
    ref. https://www.maplesoft.com/support/help/maple/view.aspx?path=MathApps%2FProjectionOfVectorOntoPlane'''

    nn = normv(n)
    l_x = normv(np.cross(y_dir, nn))
    l_y = normv(np.cross(l_x, nn))
    u_proj = u - np.dot(u,nn) / np.linalg.norm(nn)**2 * nn

    # print('proj x',)
    # print('proj y',)

    # u_proj /= np.linalg.norm(u_proj)
    return normv(np.hstack((-np.dot(u_proj, l_x), -np.dot(u_proj, l_y))))
    #_p2-u_proj_p1

def transform(spar_data0, add_pedestal_top=False, add_pedestal_bottom=False):
    ''' resultant data collection:\n
    angle, distance, slope vector
    '''
    v_ref = np.array([1,0,0])

    spar_data1 = np.roll(spar_data0, -1, axis=0)
    segment_vect = spar_data1 - spar_data0

    rot_P1 = np.zeros_like(spar_data0)
    rot_P2 = np.ones_like(spar_data0) * [0,0,1]
    buff_1 = []
    buff_2PS = []
    buff_2v = []
    buff_2R = []
    buff_2ang = []
    buff_slope_2v =[]
    for (sP1, sP2, sQ1, sQ2) in zip(rot_P1, rot_P2, spar_data0, spar_data1):
        buff_PS = []
        buff_v = []
        buff_R = []
        buff_ang = []
        buff_slope_v = []

        for (P1, P2, Q1, Q2) in zip(sP1, sP2, sQ1, sQ2):
            PS, slope_v, R = cross_point(P1, P2, Q1, Q2)
            v=proj_vector2plane(Q2-Q1, slope_v, np.array([0, 0, 1]))*10

            buff_PS.append(PS)
            buff_R.append(R)
            buff_v.append(v)
            buff_slope_v.append(slope_v)

        buff_2PS.append(np.vstack(buff_PS))
        buff_2R.append(np.vstack(buff_R))
        buff_2v.append(np.vstack(buff_v))
        buff_slope_2v.append(np.vstack(buff_slope_v))
#
#angle calc section
#
    a=np.array(buff_2PS)+np.array(buff_slope_2v)
    b=np.roll(a,-1,axis=0)
    rot_ang=np.zeros(a.shape[:2])
#-->calculate offset for the first idx_row (along spars)
    a_first=np.flipud(a[0])
    b_first=np.roll(a_first,-1, axis=0)

    rot_first_ang=np.zeros(a_first.shape[:1])
    # print(a_first, b_first)
    for i, (var_a, var_b) in enumerate(zip(a_first,b_first)):
        rot_first_ang[i]= angle(var_a, var_b, norm=False)

    rot_first_ang=np.cumsum(np.insert(rot_first_ang[:-1],0,0))
#--<
#-->calculate offsets between section points
    for i, (var_a, var_b) in enumerate(zip(a,b)):
        # print('break')
        for j, (var2a, var2b) in enumerate(zip(var_a, var_b)):
            rot_ang[i,j]= angle(var2a, var2b)

#--<
    # print('rot ang',np.degrees(rot_ang))
    angle_arr = rot_ang
    angle_arr = np.cumsum(np.rot90(rot_ang) , axis=1)
    angle_arr+= np.vstack(rot_first_ang)
    angle_arr = np.degrees(angle_arr)
#--<

#return ang, - rotation angle
#R,  - surface radius
#PS, - attachement point to axis [0,0,z]
#v - projected slope
    # print(angle_arr)
    # print(buff_2PS)
    # ang_arr = buff_2ang
    r_arr = np.roll(np.rot90(np.array(buff_2R  )),-1, axis=1)[:,:,0]
    z_arr = np.roll(np.rot90(np.array(buff_2PS )),-1, axis=1)[:,:,2]
    v_arr = np.roll(np.rot90(np.array(buff_2v )),-1, axis=1)

    return (angle_arr,
            r_arr,
            z_arr,
            v_arr)

def add_pedestal(pos, pedestal_params, add_pedestal_top=True, add_pedestal_bottom=True):
    r = pedestal_params
    h = np.array([0]+pedestal_params[1:])
    z1 = np.cumsum(h)
    z2 = np.cumsum(h[::-1])
    n_spars=pos.shape[1]

    if add_pedestal_bottom:
        pos[:,:,2] += np.sum(h)
        p_arr_1 = np.tile(pos[0],(3,1,1))
        p_arr_1[:,:,2] = np.vstack(z1)
        p_arr_1[0,:,0] = r
        pos=np.vstack([p_arr_1, pos])

    if add_pedestal_top:
        p_arr_2 = np.tile(pos[-1],(3,1,1))
        p_arr_2[:,:,2] += np.vstack(z1)
        p_arr_2[-1,:,0] = r
        pos=np.vstack([pos, p_arr_2])
    return pos

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
        idx=np.where(c<2)[0]
        if np.size(idx):
            # print(p)
            # print(c[idx])
            # print(np.size(p[:,0])-1)
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
    # print(p_arr.shape[0])
    # print('section: ',i)
    if p_arr.shape[0]>3:
        prof=make_loop(p_arr)
    return chain_list

def extract_global_params(dxf, layer):
    global_dict=collections.defaultdict(dict)
    for shape in dxf.entities:
        if layer in shape.layer and shape.dxftype == 'MTEXT':
            for block in shape.raw_text.split('[')[1:]:

                d_key =  re.findall('.*(\w+)\].*',block)
                if d_key:
                    key = d_key[0]
                # print('dkeu')
                # print(block)
                # print(key)
                    d_feed =  re.findall('feed *= *([.0-9]+)',block)
                    d_power = re.findall('power *= *([.0-9]+)',block)
                    d_angle = re.findall('angle *= *([-.0-9]+)',block)
                    d_radius =re.findall('radius *= *([-.0-9]+)',block)
                    d_cut_dir =re.findall('cut_dir.*=.*(c?cw).*',block)
                    d_split =re.findall('split.*=.*([0-9]+).*',block)

                    if d_feed:    global_dict[key]['feed']     = np.float(d_feed[0])
                    if d_power:   global_dict[key]['power']    = np.float(d_power[0])
                    if d_angle:   global_dict[key]['angle']    = np.float(d_angle[0])
                    if d_split:   global_dict[key]['split']    = np.int(d_split[0])
                    if d_radius:  global_dict[key]['radius']   = np.array([0,0,np.float(d_radius[0])])
                    if d_cut_dir: global_dict[key]['cut_dir']  = d_cut_dir

    return global_dict

def extract_params(dxf, layer):
    glob_prop_dict=collections.defaultdict(dict)
    forced_key=None

    for shape in dxf.entities:

        if layer in shape.layer and shape.dxftype == 'MTEXT':

            glob_prop_dict={'feed':200,
                            'ref_coord':np.zeros(3),
                            'power':0,
                            'angle':0,
                            'radius':np.zeros(3),
                            'layer':layer,
                            'split':None}

            porp_str = shape.raw_text

            d_feed =  re.findall('feed *= *([.0-9]+)', porp_str)
            d_power = re.findall('power *= *([.0-9]+)', porp_str)
            d_angle = re.findall('angle *= *([-.0-9]+)', porp_str)
            d_radius =re.findall('radius *= *([-.0-9]+)', porp_str)
            d_cut_dir =re.findall('cut_dir.*=.*(c?cw).*', porp_str)
            d_split =re.findall('split.*=.*([0-9]+).*', porp_str)
            d_coord =re.findall('.*coord_0.*', porp_str)

            glob_prop_dict['layer']=layer
            d_key =  re.findall('.*\[(\w+)\].*', shape.raw_text)
            if d_feed:  glob_prop_dict['feed']     = np.float(d_feed[0])
            if d_power: glob_prop_dict['power']    = np.float(d_power[0])
            if d_angle: glob_prop_dict['angle']    = np.float(d_angle[0])
            if d_split: glob_prop_dict['split']    = np.int(d_split[0])
            if d_radius:glob_prop_dict['radius'   ]= np.array([0,0,np.float(d_radius[0])])
            if d_cut_dir: glob_prop_dict['cut_dir']= d_cut_dir
            if d_coord: glob_prop_dict['ref_coord'] = [x for x in shape.insert]

            if d_key:
                forced_key = d_key[0]

    return glob_prop_dict, forced_key

def dxf_read_1(dxf, layer_name, dec_acc, n_arc, l_arc):
    tol = dec_acc
    knots_list = []
    elements_list = []
    hrd_knots_list=[]
    hrd_element_list=[]
    segment_bounds=[]
    start_coord=()
    line_count = 0
    arc_count = 0
    circle_count = 0
    path_offset =[0,0,0]
    cut_dir_marker=np.array([])
    cut_dir = 'ccw'
    split = 0
    struct_data_list=[]
    dxf_data = np.dtype([('segm', np.float, [2,3]),('props', np.int)])
    prop_list=[]
    elements_list = []
    glob_prop_dict= {}
    prop_dict = collections.defaultdict(dict)

    global_params = extract_global_params(dxf, 'glob_props')

    for p, layer in enumerate(layer_name):

        local_params, forced_key = extract_params(dxf, layer)

        if local_params:
            prop_dict[p] = local_params
        else:
            default_params={'feed':200,
            'ref_coord':np.zeros(3),
            'power':0,
            'angle':0,
            'radius':np.zeros(3),
            'layer':layer,
            'split':None}

            prop_dict[p] = default_params

        for shape in dxf.entities:

            if shape.layer == layer:
                if shape.dxftype == 'MTEXT':
                    if 'start'in shape.raw_text:
                        # print('found start')
                        start_coord = tuple(round(x, tol) for x in shape.insert)

                if shape.dxftype == 'LINE':
                    line_count += 1
                    p1 = tuple(round(x, tol) for x in shape.start)
                    p2 = tuple(round(x, tol) for x in shape.end)
                    if p1!=p2:
                        knots_list.append(p1)
                        knots_list.append(p2)
                        elements_list.append([p1, p2])
                        prop_list.append(p)

                if shape.dxftype == 'ARC':
                    arc_count += 1
                    ARC_knots_list = []
                    n = n_arc  # number of segments
                    min_len = l_arc
                    O = shape.center
                    R = shape.radius
                    angl_1 = shape.start_angle * np.pi / 180
                    angl_2 = shape.end_angle * np.pi / 180

                    if angl_2 >= angl_1:
                        angl_list = np.linspace(angl_1, angl_2, n)
                    else:
                        angl_list = np.linspace(angl_1, angl_2 + 2 * np.pi, n)

                    arc_len = R * np.absolute(angl_2 - angl_1)

                    if arc_len / n < min_len:
                        n = max(int(arc_len / min_len), 3)

                    for angl in angl_list:
                        ARC_knots_list.append(
                            (round(O[0] + R * np.cos(angl), tol), round(O[1] + R * np.sin(angl), tol), O[2]))

                    for i in range(n - 1):
                        elements_list.append(ARC_knots_list[i:i + 2])
                        prop_list.append(p)

    element_arr=np.array(elements_list)
    prop_arr=np.array(prop_list)
    start_coord_arr = np.array(start_coord)

    return (element_arr, prop_arr, prop_dict, (start_coord_arr, cut_dir_marker, cut_dir, split))

def dxf_read_2(dxf_name, lname):

    seg_dt = np.dtype([('seg', float, (2,3)), ('prp', int, 1)])

    dwg = ezdxf.readfile(dxf_name)
    modelspace = dwg.modelspace()

    lines = modelspace.query('LINE[layer=="{}"]'.format(lname))
    arcs = modelspace.query('ARC[layer=="{}"]'.format(lname))
    mtext = modelspace.query('MTEXT[layer=="{}"]'.format(lname))

    seg_lin = np.stack([np.array([ ( np.round(( line.dxf.start, line.dxf.end ),4),hash(lname))], dtype = seg_dt) for line in lines], axis = 0)
    # print(seg_lin)

    return seg_lin
    # print(np.vstack(seg_lin))
    # for p, layer in enumerate(layer_name):
    #
    #     # local_params, forced_key = extract_params(dxf, layer)
    #
    #     # if local_params:
    #     #     prop_dict[p] = local_params
    #     # else:
    #     #     default_params={'feed':200,
    #     #     'ref_coord':np.zeros(3),
    #     #     'power':0,
    #     #     'angle':0,
    #     #     'radius':np.zeros(3),
    #     #     'layer':layer,
    #     #     'split':None}
    #     #
    #     #     prop_dict[p] = default_params
    #
    #     for shape in dxf.entities:
    #
    #         if shape.layer == layer:
    #
    #             # if shape.dxftype == 'MTEXT':
    #             #     if 'start'in shape.raw_text:
    #             #         # print('found start')
    #             #         start_coord = tuple(round(x, tol) for x in shape.insert)
    #
    #             if shape.dxftype == 'LINE':
    #                 line_count += 1
    #                 p1 = tuple(round(x, tol) for x in shape.start)
    #                 p2 = tuple(round(x, tol) for x in shape.end)
    #                 if p1!=p2:
    #                     knots_list.append(p1)
    #                     knots_list.append(p2)
    #                     elements_list.append([p1, p2])
    #                     prop_list.append(p)
    #
    #                     if shape.dxftype == 'ARC':
    #                         arc_count += 1
    #                         ARC_knots_list = []
    #                         n = n_arc  # number of segments
    #                         min_len = l_arc
    #                         O = shape.center
    #                         R = shape.radius
    #                         angl_1 = shape.start_angle * np.pi / 180
    #                         angl_2 = shape.end_angle * np.pi / 180
    #
    #                         if angl_2 >= angl_1:
    #                             angl_list = np.linspace(angl_1, angl_2, n)
    #                         else:
    #                             angl_list = np.linspace(angl_1, angl_2 + 2 * np.pi, n)
    #
    #                             arc_len = R * np.absolute(angl_2 - angl_1)
    #
    #                             if arc_len / n < min_len:
    #                                 n = max(int(arc_len / min_len), 3)
    #
    #                                 for angl in angl_list:
    #                                     ARC_knots_list.append(
    #                                     (round(O[0] + R * np.cos(angl), tol), round(O[1] + R * np.sin(angl), tol), O[2]))
    #
    #                                     for i in range(n - 1):
    #                                         elements_list.append(ARC_knots_list[i:i + 2])
    #                                         prop_list.append(p)
    #
    #                                         element_arr=np.array(elements_list)
    #                                         prop_arr=np.array(prop_list)
    #                                         start_coord_arr = np.array(start_coord)
    # return (element_arr, prop_arr, prop_dict, (start_coord_arr, cut_dir_marker, cut_dir, split))



def find_segment1(data, pt):
    sol_idx = np.array([])
    A=np.product(data[:,0,:]==pt, axis =1)
    B=np.product(data[:,1,:]==pt, axis =1)

    Aidx = np.where(A==1)[0]
    Bidx = np.where(B==1)[0]

    if Aidx.size:
        sol_idx = np.append(sol_idx, Aidx)
    if Bidx.size:
        sol_idx = np.append(sol_idx, Bidx)

    return sol_idx.astype(int)

def find_unique_segments(data,pt):
    sol_idx = np.array([])
    d_list =[]
    for line in data:
        d, i = find_nearest(line, pt)
        if i:
            d_list.append(np.dot(*line[::-1]))
        else:
            d_list.append(np.dot(*line))
    print(d_list)

def data_fix(data):
    # for line in data:
    #     print(line)
    print('data to fix', data.shape)
    d=data.reshape(-1,3)
    # for line in d:
    #     print(line)
    val, counts = np.unique(d, axis=0, return_counts=True)
    print('single',np.where(counts == 1)[0])
    idx3 = np.where(counts >=3)[0]
    print('triple and more',idx3)
    # print('counts',)
    for i, idx in enumerate(idx3):
        print(i,' #find segm    ents with:', val[idx])
        seg_idx_arr = find_segment1(data,val[idx])
        print(data[seg_idx_arr])
        print(find_unique_segments(data[seg_idx_arr], val[idx]))
    # print('unique points',np.unique(d, axis=0).shape)
    # print(data[-1])
    return 0

def extract_dxf_layer_data(fname_dxf, layer, dec_acc=4, n_arc=10, l_arc=0.1):
    dxf = dxfgrabber.readfile(fname_dxf, {"assure_3d_coords": True})
    geom, prop, prop_dict, glob_params = dxf_read_1(dxf, [layer], dec_acc, n_arc, l_arc)
    return {'geom':geom, 'prop':prop, 'prop_dict':prop_dict, 'glob_params':glob_params}


def extract_dxf_path(fname_dxf, layer_list, dec_acc=4, n_arc=10, l_arc=0.1):
    dxf = dxfgrabber.readfile(fname_dxf, {"assure_3d_coords": True})
    struct_data, prop_data, prop_dict, glob_params = dxf_read_1(dxf, layer_list, dec_acc, n_arc, l_arc)
    start_coord_arr, cut_dir_marker, cut_dir, split = glob_params
    io_path, io_rest, io_path_prop, io_rest_prop = find_io_path(struct_data, prop_data, start_coord_arr, prop_dict = prop_dict)
    pt0 = io_path[-1,-1]

    if io_rest:
        lo_path, lo_rest, lo_path_prop, lo_rest_prop = find_lo_path(io_rest, io_rest_prop, pt0, cut_dir = cut_dir, cut_dir_marker = cut_dir_marker, prop_dict = prop_dict)
    else:
        lo_path = np.array([])
        lo_rest =  np.array([])
        lo_path_prop =np.array([])
        lo_rest_prop =np.array([])

    return io_path, lo_path, io_path_prop, lo_path_prop, prop_dict

def find_nearest(arr, pt0):
    tree = spatial.cKDTree(arr)
    d, i = tree.query(pt0, k=[1])
    return d, i

def transform_pt(p, ref_coord = np.zeros(3), r=np.zeros(3), angle=0):
    # def M(axis, theta):
    #     return expm(np.cross(np.eye(3), axis/norm(axis)*np.radians(angle)))
    def M(axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
        """
        axis = np.asarray(axis)
        axis = axis/np.sqrt(np.dot(axis, axis))
        a = np.cos(theta/2.0)
        b, c, d = -axis*np.sin(theta/2.0)
        aa, bb, cc, dd = a*a, b*b, c*c, d*d
        bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
        return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                         [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                         [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

    axis = np.array([0,0.1,0])
    M0 = M(axis, np.radians(angle))
    p_abs = p+r-ref_coord
    return np.dot(M0, p_abs) #+ np.dot(M0, r)
    #+ref_coord

def plot_path(path_list, prop_list, pd_list):
    color_list = ['red', 'green', 'blue', 'yellow']
    # print(path_list)
    path_arr = np.array(path_list)
    prop_arr = np.array(prop_list)

    if len(prop_arr.shape) == 1:
        ax = plt.subplot(1,1,1)
        for j in np.arange(len(path_list)):
            for i in np.arange(len(path_list[j])):
                # print(color_list[prop_list[j][i]])
                ax.plot(path_list[j][i][:,0],path_list[j][i][:,1], color=color_list[prop_list[j][i]])
        #
        plt.grid(True)
        plt.show()
def make_path_chain(io,lo):
    return np.vstack((io[:,0,:],    io[-1,1,:],
                      lo[:,0,:],    lo[-1,1,:],
                      io[::-1,1,:], io[0,0,:]))

def make_prop_chain(io, lo):
    return np.hstack((io, io[-1], lo, lo[-1]))


def make_path_chain_single(io):
    return np.vstack((io[:,0,:],    io[-1,1,:]))

def make_prop_chain_single(io):
    return np.hstack((io, io[-1]))

def make_projection(chain1, chain2, chain1_prop, chain2_prop, prop_dict1, prop_dict2, ret_full = False):
    offset_list1 = []
    offset_list2 = []
    radius_list1 = []
    radius_list2 = []
    # print(prop_dict1)
    for pr1, pr2 in zip(chain1_prop, chain2_prop):

        # print(prop_dict1[pr1]['ref_coord'])
        offset_list1.append(prop_dict1[pr1]['ref_coord'])
        offset_list2.append(prop_dict2[pr2]['ref_coord'])
        radius_list1.append(prop_dict1[pr1]['radius'])
        radius_list2.append(prop_dict2[pr2]['radius'])

    # print(chain2_prop)
    # print(np.array(offset_list1))
    ch1 = chain1 - np.array(offset_list1) + np.array(radius_list1)
    ch2 = chain2 - np.array(offset_list2) + np.array(radius_list2)

    proj_ch1 = []
    proj_ch2 = []
    for var1, var2 in zip(ch1, ch2):
        # print(var2-var1)
        proj_ch1.append(line_plane_intersect(var2, var1, np.array([0,0,-310]), np.array([0,0,1])))
        proj_ch2.append(line_plane_intersect(var2, var1, np.array([0,0, 310]), np.array([0,0,1])))
        # proj_ch2 = cncfclib.proj_vector2plane(, n, y_dir = np.array([0,0,1])):

    # print('ch1',ch1[-1])
    proj_ch1_arr = np.array(proj_ch1)
    # print('proj',proj_ch1_arr[-1])


    # print('ch2',ch2[-1])
    proj_ch2_arr = np.array(proj_ch2)
    # print('proj',proj_ch2_arr[-1])

    if ret_full:
        return proj_ch1_arr, proj_ch2_arr
    else:
        return proj_ch1_arr[:,:2], proj_ch2_arr[:,:2]

def transform_chain(chain, pd1):
    chain_tr = []
    for pt in chain:
        chain_tr.append(transform_pt(pt,ref_coord = pd1['ref_coord'], r=pd1['radius'], angle=pd1['angle']))
    return np.array(chain_tr)

def plot_path1(ss):
    color_list = ['red', 'green', 'blue', 'yellow']
    # print(path_list)
    # ax = plt.subplot(1,1,1)
    fig = plt.figure()
    ax = fig.gca(projection='3d')

### quick fix to get projections
###
###
    for i, var1 in enumerate(ss):
        # print('VAR ',var1)
        if var1:
            io_path1, lo_path1, io_prop1, lo_prop1, prop_dict1 = var1[0]
            io_path2, lo_path2, io_prop2, lo_prop2, prop_dict2 = var1[1]

            # gcode.append('{0[0]}layers: {1[0]} {1[1]}{0[1]}'.format(gcode_conf['comment'],[prop_dict1[0]['layer'], prop_dict2[0]['layer']]))
            # gcode.append('{0[0]}sequence :{1}{0[1]}'.format(gcode_conf['comment'],i))

            if lo_path1 and lo_path2:
                chain1_path = make_path_chain(io_path1, lo_path1)
                chain2_path = make_path_chain(io_path2, lo_path2)
                chain1_prop = make_prop_chain(io_prop1, lo_prop1)
                chain2_prop = make_prop_chain(io_prop2, lo_prop2)
            else:
                chain1_path = make_path_chain_single(io_path1)
                chain2_path = make_path_chain_single(io_path2)
                chain1_prop = make_prop_chain_single(io_prop1)
                chain2_prop = make_prop_chain_single(io_prop2)

            p1tproj, p2tproj = make_projection(chain1_path, chain2_path, chain1_prop, chain2_prop, prop_dict1, prop_dict2, ret_full = True)
            pd1 = prop_dict1[0]
            pd2 = prop_dict2[0]
            print(pd1)
            p1tproj = transform_chain(p1tproj, pd1)#ref_coord = pd1['ref_coord'], r=pd1['radius'], angle=pd1['angle'])
            p2tproj = transform_chain(p2tproj, pd2)#ref_coord = pd2['ref_coord'], r=pd2['radius'], angle=pd2['angle'])
            ax.plot(p1tproj[:,0],p1tproj[:,1],p1tproj[:,2], color='k')
            ax.plot(p2tproj[:,0],p2tproj[:,1],p2tproj[:,2], color='k')
###
###
###





    for seq in ss:
        # print(seq)
        for pp in seq:
            # print('TTTTTTTTTTTTTTt')
            # print(pp)
            # for plane_pp in pp:
            io_path1, lo_path1, io_path_prop1, lo_path_prop1, prop_dict1 = pp
            for path, prop in [(io_path1, io_path_prop1), (lo_path1, lo_path_prop1)]:

                for i in np.arange(path.shape[0]):
                    p1 = path[i][0]
                    p2 = path[i][1]

                    # p1_proj, p2_proj = cncfclib.make_projection(chain1_path, chain2_path, chain1_prop, chain2_prop, prop_dict1, prop_dict2)

                    prop_n = prop[i]

                    pd = prop_dict1[prop_n]
                    c = color_list[prop_n]

                    p1t = transform_pt(p1,ref_coord = pd['ref_coord'], r=pd['radius'], angle=pd['angle'])
                    p2t = transform_pt(p2,ref_coord = pd['ref_coord'], r=pd['radius'], angle=pd['angle'])
                    # if p1t[2] < 0  or p2t[2] <0:
                        # print(p1t[2], p2t[2])
                    # px = np.hstack((p1t[0], p2t[0]))
                    # py = np.hstack((p1t[1], p2t[1]))
                    # pz = np.hstack((p1t[2], p2t[2]))
                    ax.quiver(p1t[0], p1t[1], p1t[2],
                              p2t[0]-p1t[0], p2t[1]-p1t[1], p2t[2]-p1t[2], color=c, arrow_length_ratio=0)
    plt.grid(True)
    plt.show()

def find_segments(arr, member_pt):
    idx = np.where((arr[:,0,:] == member_pt) | (arr[:,1,:] == member_pt))[0]
    # print('member point',member_pt)
    # print(arr[:,0,:])
    # print(idx)
    return arr[idx,:,:]

def sort_segments(arr, start_pt, stop_pt=np.array([]), close_loop = False, return_idx=False, prop_data=np.array([])):
    '''Make a /chain sorting/ of segments defined in the arr. The starting point is definied by start_pt. If there is no exact pt matching, the function uses the nearest end pt
    PARAMETERS:
        arr = np.array( nx2x3 ) ex: [[[x,y,z],[x1,y1,z1]],[[x2,y2,y3],[x1,y1,z1]],...]
        start_pt = np.array([x,y,z]) - vector
        stop_pt = np.array([x,y,z]) - vector, default empty.
        close_loop = False|True, default False.
        return_idx = False|True, default False. If True the function returns seorted indices of arr
    RETURNS:
        sol - sorted segments
        rest - remaining segments. in general subtract sol from arr
    '''
    rest=np.array([])
    #split segments array to 2 columns
    A_arr = arr[:,0,:]
    B_arr = arr[:,1,:]
    C_arr = prop_data
    #pt0 is a loop variable
    pt0 = start_pt

    # print(prop_data)

        #make data storage buffers
    sorted_arr = np.zeros_like(arr)
    sorted_prop_arr = np.zeros_like(prop_data)
    rest = np.array([[],[]])
    rest_prop = np.array([])

    for i in np.arange( np.shape(arr)[0]):
        #find the nearest point to pt0
        #return X_d - distance; X_i - position in the arr

        A_d, A_i = find_nearest(A_arr, pt0)
        B_d, B_i = find_nearest(B_arr, pt0)

        #If a point in A_arr is closer to pt0
        if A_d[0]<B_d[0]:

            sorted_arr[i,0,:]= A_arr[A_i[0]]
            sorted_arr[i,1,:]= B_arr[A_i[0]]
            if prop_data.size>0:
                sorted_prop_arr[i]= C_arr[A_i[0]]

            pt0 = B_arr[A_i[0]]
            A_arr = np.delete(A_arr, A_i[0], axis = 0)
            B_arr = np.delete(B_arr, A_i[0], axis = 0)
            if prop_data.size>0:
                C_arr = np.delete(C_arr, A_i[0], axis = 0)

        else:

            sorted_arr[i,0,:]= B_arr[B_i[0]]
            sorted_arr[i,1,:]= A_arr[B_i[0]]
            if prop_data.size>0:
                sorted_prop_arr[i]= C_arr[B_i[0]]

            pt0 = A_arr[B_i[0]]

            A_arr = np.delete(A_arr, B_i[0], axis = 0)
            B_arr = np.delete(B_arr, B_i[0], axis = 0)
            if prop_data.size>0:
                C_arr = np.delete(C_arr, B_i[0], axis = 0)

        #if stop_pt within current segment then break
        if np.array_equal(stop_pt, pt0):
            rest = np.stack((A_arr,B_arr), axis=1)
            rest_prop = C_arr
            break
# print(sorted_arr)
    #close loop solution array modyfication
    if prop_data.size>0:
        if close_loop:
            sol_io = np.vstack((sorted_arr[:i+1],
                             np.array([ [sorted_arr[-1,1,:], sorted_arr[0,0,:]] ])))
            sol_prop = np.vstack((sorted_prop_arr[:i+1], sorted_prop_arr[-1]))
        else:
            sol_io = sorted_arr[:i+1]
            sol_prop = sorted_prop_arr[:i+1]

        sol = (sol_io, rest, sol_prop, rest_prop)

    else:
        if close_loop:
            sol_io = np.vstack((sorted_arr[:i+1],
                            np.array([ [sorted_arr[-1,1,:], sorted_arr[0,0,:]] ])))

        else:
            sol_io = sorted_arr[:i+1]

        sol = (sol_io, rest)

    return sol

def sort_segments1(seg_arr, start_pt, stop_pt=np.array([]), close_loop = False, return_idx=False, prop_data=np.array([])):
    '''Make a /chain sorting/ of segments defined in the arr. The starting point is definied by start_pt. If there is no exact pt matching, the function uses the nearest end pt
    PARAMETERS:
        arr = np.array( nx2x3 ) ex: [[[x,y,z],[x1,y1,z1]],[[x2,y2,y3],[x1,y1,z1]],...]
        start_pt = np.array([x,y,z]) - vector
        stop_pt = np.array([x,y,z]) - vector, default empty.
        close_loop = False|True, default False.
        return_idx = False|True, default False. If True the function returns seorted indices of arr
    RETURNS:
        sol - sorted segments
        rest - remaining segments. in general subtract sol from arr
    '''
    def FlipSeg(seg_arr):
        seg_arr['seg']=seg_arr['seg'][::-1]
        return seg_arr

    rest=np.array([])
    #split segments array to 2 columns
    A_arr = seg_arr['seg'][:,0,:]
    B_arr = seg_arr['seg'][:,1,:]
    seg_arr_rest = deepcopy(seg_arr)
    # C_arr = prop_data
    #pt0 is a loop variable
    pt0 = start_pt

    # print(prop_data)

    #make data storage buffers
    sorted_arr = np.zeros_like(seg_arr)
    # rest = np.array([[],[]])
    # rest_prop = np.array([])

    for i in np.arange( np.shape(seg_arr)[0]):
        #find the nearest point to pt0
        #return X_d - distance; X_i - position in the arr

        A_d, A_i = find_nearest(A_arr, pt0)
        B_d, B_i = find_nearest(B_arr, pt0)

        #If a point in A_arr is closer to pt0
        if A_d[0]<B_d[0]:
            sorted_arr[i] = seg_arr[A_i[0]]
            pt0 = B_arr[A_i[0]]
            A_arr = np.delete(A_arr, A_i[0], axis = 0)
            B_arr = np.delete(B_arr, A_i[0], axis = 0)
        else:
            sorted_arr[i]= FlipSeg(seg_arr[B_i[0]])
            pt0 = A_arr[B_i[0]]
            A_arr = np.delete(A_arr, B_i[0], axis = 0)
            B_arr = np.delete(B_arr, B_i[0], axis = 0)

        #if stop_pt within current segment then break
        if np.array_equal(stop_pt, pt0):
            rest = np.stack((A_arr,B_arr), axis=1)
            rest_prop = C_arr
            break
    #close loop solution array modyfication
    if close_loop:
        sol_io = np.vstack((sorted_arr[:i+1],
                        np.array([ [sorted_arr[-1,1,:], sorted_arr[0,0,:]] ])))

    else:
        sol_io = sorted_arr[:i+1]

    sol = (sol_io, rest)

    return sol



def sort_loop(arr, start_pt, dir='cw'):
    idx = find3dpoint(arr, start_pt)
    # v1
    # v2
    sol, rest = sort_segments(arr, start_pt, close_loop = True)
    return sol, rest

def make_split(arr, props, splits):
    v = arr[:,1] - arr[:,0]

    l_norm = np.linalg.norm(v, axis=1)
    l_cs = np.hstack((0, np.cumsum(l_norm)))

    n_seg = np.linspace(0,np.sum(l_norm),splits)
    x=np.interp(n_seg, l_cs, np.hstack((arr[:,0,0], arr[-1,1,0])))
    y=np.interp(n_seg, l_cs, np.hstack((arr[:,0,1], arr[-1,1,1])))
    z=np.interp(n_seg, l_cs, np.hstack((arr[:,0,2], arr[-1,1,2])))

    arr_buff0 = np.column_stack((x, y, z))
    arr_buff1= np.roll(arr_buff0, -1, axis=0)
    arr_buff  = np.stack((arr_buff0, arr_buff1),axis=1)[:-1]
    prop_buff = np.ones(arr_buff.shape[0], dtype=np.int) * props[0]

    return arr_buff, prop_buff

def MakeSplit(seg_arr, n, args={}):
    #seg_arr - segments array
    #n - number of new segments
    v = arr[:,1] - arr[:,0]

    l_norm = np.linalg.norm(v, axis=1)
    l_cs = np.hstack((0, np.cumsum(l_norm)))

    n_seg = np.linspace(0,np.sum(l_norm),splits)
    x=np.interp(n_seg, l_cs, np.hstack((arr[:,0,0], arr[-1,1,0])))
    y=np.interp(n_seg, l_cs, np.hstack((arr[:,0,1], arr[-1,1,1])))
    z=np.interp(n_seg, l_cs, np.hstack((arr[:,0,2], arr[-1,1,2])))

    arr_buff0 = np.column_stack((x, y, z))
    arr_buff1= np.roll(arr_buff0, -1, axis=0)
    arr_buff  = np.stack((arr_buff0, arr_buff1),axis=1)[:-1]
    prop_buff = np.ones(arr_buff.shape[0], dtype=np.int) * props[0]

    return arr_buff, prop_buff

def find_io_path(arr, prop_data, start_pt=np.array([]), return_idx=False, prop_dict={}):
    '''finds open loop chain. the stop point is at T-junction or the free end
        PARAMETERS:
            arr = np.array( nx2x3 ) ex: [[[x,y,z],[x1,y1,z1]],[[x2,y2,y3],[x1,y1,z1]],...]
            start_pt = np.array([x,y,z]) - vector
        RETURNS:
            sol - sorted segments
            rest - remaining segments. in general subtract sol from arr
    '''

    knots_arr = arr.reshape(-1,3)

    unique_knots, counts = np.unique(knots_arr, return_counts=True, axis=0)

    unique_knots_1 = unique_knots[np.where(counts == 1)[0]]
    unique_knots_3 = unique_knots[np.where(counts >= 3)[0]]

    if start_pt.size:
        IO_knot_d, IO_knot_i = find_nearest(unique_knots_1, start_pt)
        IO_knot = unique_knots_1[IO_knot_i]
        stop_knot = np.array([])
    else:
        IO_knot = unique_knots_1
        stop_knot = unique_knots_3[0]

    sol, rest, sol_prop, rest_prop = sort_segments(arr, IO_knot[0], stop_pt=stop_knot, prop_data = prop_data)

    if prop_dict:
        for p in set(sol_prop):
            if prop_dict[p]['split']:
                sol_list = []
                prop_list = []
                splits = prop_dict[p]['split']
                p_idxs = np.where(sol_prop == p)[0]
                q_idxs = np.where(sol_prop != p)[0]
                # print('ZZZZZZZZZZZZZZZZZz', p_idxs)
                continous_ranges_list = []
                ranges_to_split_list = []
                ranges_to_save_list = []

                for k, g in groupby(enumerate(p_idxs.tolist()), lambda ix: ix[0]-ix[1]):
                    ranges_to_split = list(map(itemgetter(1), g))
                    ranges_to_split_list.append(ranges_to_split)

                for k, b in groupby(enumerate(q_idxs.tolist()), lambda ix: ix[0]-ix[1]):
                    ranges_to_save = list(map(itemgetter(1), b))
                    ranges_to_save_list.append(ranges_to_save)

                A = ranges_to_split_list
                B = ranges_to_save_list

                if not B:
                    # print('no ranges to save')
                    for split_group in A:
                        split_buff, split_prop_buf = make_split(sol[split_group], sol_prop[split_group], splits)
                    sol_list.append(split_buff)
                    prop_list.append(split_prop_buf)
                else:

                    if len(A)==len(B):

                        for split_group, save_group in zip(A, B):
                            split_buff, split_prop_buf = make_split(sol[split_group], sol_prop[split_group], splits)
                            if A[0][0]<B[0][0]:
                                sol_list.append(split_buff)
                                sol_list.append(sol[save_group])
                                prop_list.append(split_prop_buf)
                                prop_list.append(sol_prop[save_group])
                            else:
                                sol_list.append(sol[save_group])
                                sol_list.append(split_buff)
                                prop_list.append(sol_prop[save_group])
                                prop_list.append(split_prop_buf)

                    else:
                        if len(A)>len(B):
                            for split_group, save_group in zip(A[:-1],B):
                                split_buff, split_prop_buf = make_split(sol[split_group], sol_prop[split_group], splits)
                                sol_list.append(split_buff)
                                sol_list.append(sol[save_group])
                                prop_list.append(split_prop_buf)
                                prop_list.append(sol_prop[save_group])
                            split_buff, split_prop_buf = make_split(sol[A[-1]], sol_prop[A[-1]], splits)
                            sol_list.append(split_buff)
                            prop_list.append(split_prop_buf)

                        else:
                            for split_group, save_group in zip(A,B[:-1]):
                                split_buff, split_prop_buf = make_split(sol[split_group], sol_prop[split_group], splits)
                                sol_list.append(sol[save_group])
                                sol_list.append(split_buff)
                                prop_list.append(sol_prop[save_group])
                                prop_list.append(split_prop_buf)
                            sol_list.append(sol[B[-1]])
                            prop_list.append(sol_prop[B[-1]])

                sol = np.concatenate(sol_list, axis=0)
                sol_prop = np.hstack(prop_list)

    return sol, rest, sol_prop, rest_prop

def find_io_path1(arr, start_pt=np.zeros(3), prp_dict={}, return_idx=False, prop_dict={}):
    '''finds open loop chain. the stop point is at T-junction or the free end
        PARAMETERS:
            arr = np.array( nx2x3 ) ex: [[[x,y,z],[x1,y1,z1]],[[x2,y2,y3],[x1,y1,z1]],...]
            start_pt = np.array([x,y,z]) - vector
        RETURNS:
            sol - sorted segments
            rest - remaining segments. in general subtract sol from arr
    '''
    print(prp_dict)
    if prp_dict:
        if 'start' in prp_dict.keys():

            start = prp_dict[arr['prp'][0]]['start']

    knots_arr = arr['seg'].reshape(-1,3)
    # print(knots_arr)
    unique_knots, counts = np.unique(knots_arr, return_counts=True, axis=0)

    unique_knots_1 = unique_knots[np.where(counts == 1)[0]]
    # print(unique_knots_1)
    unique_knots_3 = unique_knots[np.where(counts >= 3)[0]]
    # print('u knots 3',unique_knots_3)

    if start_pt.size:
        IO_knot_d, IO_knot_i = find_nearest(unique_knots_1, start_pt)
        IO_knot = unique_knots_1[IO_knot_i]
        stop_knot = np.array([])
    else:
        IO_knot = unique_knots_1
        stop_knot = unique_knots_3[0]

    sol, rest = sort_segments1(arr, start)

    print(sol)

    return sol, rest




def find_lo_path(arr, prop_data, start_pt, return_idx=False, close_loop = True, cut_dir = 'ccw', cut_dir_marker = np.array([]), prop_dict={}):
    dir_dict={'ccw':1, 'cw':-1}
    # print(arr)
    # print(start_pt)

    z = find3dpoint(arr, start_pt)
    # print('z:', z)
    s1=np.array([arr[z[0]]])
    s2=np.array([arr[z[1]]])
    sp=start_pt

    s1_path, s1_rest = sort_segments(s1, sp)
    s2_path, s2_rest = sort_segments(s2, sp)

    if cut_dir_marker.size:
        if np.linalg.norm(cut_dir_marker-s1_path[-1])<=np.linalg.norm(cut_dir_marker-s2_path[-1]):
            lo_rest=np.delete(arr,z[1], axis=0)
            lo_rest_prop=np.delete(prop_data,z[1])
            removed_seg = arr[z[1]]
            removed_prop = prop_data[z[1]]
        else:
            lo_rest=np.delete(arr,z[0], axis=0)
            lo_rest_prop=np.delete(prop_data,z[0])
            removed_seg = arr[z[0]]
            removed_prop = prop_data[z[0]]
    else:
        u=np.diff(s1_path[0], axis=0)[0]
        v=np.diff(s2_path[0], axis=0)[0]
        loop_dir = np.cross(u,v)[2]
        # print(cut_dir)
        if dir_dict[cut_dir]*loop_dir>=0:
            lo_rest=np.delete(arr,z[1], axis=0)
            lo_rest_prop=np.delete(prop_data,z[1])
            removed_seg = arr[z[1]]
            removed_prop = prop_data[z[1]]
        else:
            lo_rest=np.delete(arr,z[0], axis=0)
            lo_rest_prop=np.delete(prop_data,z[0])
            removed_seg = arr[z[0]]
            removed_prop = prop_data[z[0]]

    # print('input arr',arr, arr.shape[0])
    # print('sp',sp)
    sol, rest, sol_prop, rest_prop = sort_segments(lo_rest, sp, prop_data = lo_rest_prop)
    if close_loop:
        # print('TTTTTTTTTTTTt')
        # print(sol[-1][1])
        # print(removed_seg)

        removed_seg_sorted, dummy = sort_segments(np.array([removed_seg]), sol[-1][1])
        # print(removed_seg_sorted)
        sol = np.vstack((sol, removed_seg_sorted))
        # print(sol_prop, np.array([removed_prop]))
        sol_prop = np.hstack((sol_prop, np.array([removed_prop])))
    return sol, rest, sol_prop, rest_prop

#DXF2KNOTS functions potentialy to remove
def sort2match(arr, p):
    if np.array_equal(s1[0], sp):
        v=np.diff(s1, axis=0)
        # print('norm')
    else:
        v=np.diff(s1[::-1])
        # print('rev')
    return v

def sub_points(p1, p2):
    vect = []
    p1 = [x for x in p1[0]]
    p2 = [x for x in p2[0]]

    if len(p1) == len(p2):
        for i, n in enumerate(p2):
            vect.append(n - p1[i])
        return vect
    return len(p1) * [None]


def knots_rank_find(knots_rank, rank):
    knots = [x[0] for x in knots_rank if x[1] == rank]
    if len(knots) > 0:
        return knots
    else:
        return [None]


def knots_rank_list(el_kt_list, sorted_knots, skip_knot):
    knots_rank = []
    for var in sorted_knots:
        if var[0] == skip_knot:
            knots_rank.append([var[0], None])
        else:
            knots_rank.append([var[0], [x for a in el_kt_list for x in a].count(var[0])])
    return knots_rank


def knot2coord(sorted_knots, knot):
    for var in sorted_knots:
        if var[0] == knot:
            return var[1]
    return None


def knots2file(name, io_path, sorted_knots):
    f = open(name, 'w')

    for var in io_path:
        coord = knot2coord(sorted_knots, var[0])
        f.write('{0:.3f} {1:.3f}\n'.format(coord[0], coord[1]))

    coord = knot2coord(sorted_knots, io_path[-1][1])
    f.write('{0:.3f} {1:.3f}\n'.format(coord[0], coord[1]))

    f.close()

def ct_len(io_path, sorted_knots):
    coord_list = []
    for var in io_path:
        coord = knot2coord(sorted_knots, var[0])
        coord_list.append((coord[0], coord[1]))

    coord = knot2coord(sorted_knots, io_path[-1][1])
    coord_list.append((coord[0], coord[1]))

    coord_arr=np.array(coord_list)
    l_arr = np.linalg.norm(np.diff(coord_arr,axis=0), axis=1)
    return np.sum(l_arr)

def knots2file_1(name, section_list):
    f = open(name, 'w')
    # print('TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT')
    for var in section_list:
        # coord = knot2coord(sorted_knots, var[0])
        # print(var)
        f.write('{0:.3f} {1:.3f} {2:.3f}\n'.format(var[0], var[1], z_coord))

    # coord = knot2coord(sorted_knots, io_path[-1][1])
    # f.write('{0:.3f} {1:.3f}\n'.format(coord[0], coord[1]))

    f.close()


def knots_dict(knots_list):
    return [[i, var] for i, var in enumerate(list(set(knots_list)))]

def elements_coords2knots(el_list, kt_list):
    el_kt_list = []
    for el in el_list:
        for kt in kt_list:
            if kt[1] == el[0]:
                p1 = kt[0]
            if kt[1] == el[1]:
                p2 = kt[0]
        el_kt_list.append([p1, p2])
    return el_kt_list


def elements_knots2coords(el_list, kt_list):
    el_coord_list = []
    for el in el_list:
        for kt in kt_list:
            if kt[0] == el[0]:
                p1 = kt[1]
            if kt[0] == el[1]:
                p2 = kt[1]
        el_coord_list.append([p1, p2])
    return el_coord_list


def find_path(crit, el_kt_list, sorted_knots, excl_knot):
    path = []

    knots_rank = knots_rank_list(el_kt_list, sorted_knots, excl_knot)

    curr_knot = knots_rank_find(knots_rank, 1)
    last_knot = knots_rank_find(knots_rank, 3)

    curr_element=[]

    # print '\nfinding path'
    while not ((curr_element is None) or curr_knot[0]==last_knot[0]):
        # print '\rpool size: {0}'.format(len(path)),

        curr_element=next((element for element in el_kt_list if curr_knot[0] in element), None)
        if not (curr_element is None):

            if curr_element[0] == curr_knot[0]:
                curr_knot=[curr_element[1]]
                path.append(curr_element)
            else:
                curr_knot=[curr_element[0]]
                path.append(curr_element[::-1])

            el_kt_list.remove(curr_element)

    if crit == 1:
        path.append([path[-1][1], path[0][0]])
    # print '\n'
    return path

def extract_params(dxf, layer):
    glob_prop_dict=collections.defaultdict(dict)
    forced_key=None

    for shape in dxf.entities:

        if layer in shape.layer and shape.dxftype == 'MTEXT':

            glob_prop_dict={'feed':200,
                            'ref_coord':np.zeros(3),
                            'power':0,
                            'angle':0,
                            'radius':np.zeros(3),
                            'layer':layer,
                            'split':None}

            porp_str = shape.raw_text

            d_feed =  re.findall('feed *= *([.0-9]+)', porp_str)
            d_power = re.findall('power *= *([.0-9]+)', porp_str)
            d_angle = re.findall('angle *= *([-.0-9]+)', porp_str)
            d_radius =re.findall('radius *= *([-.0-9]+)', porp_str)
            d_cut_dir =re.findall('cut_dir.*=.*(c?cw).*', porp_str)
            d_split =re.findall('split.*=.*([0-9]+).*', porp_str)
            d_coord =re.findall('.*coord_0.*', porp_str)

            glob_prop_dict['layer']=layer
            d_key =  re.findall('.*\[(\w+)\].*', shape.raw_text)
            if d_feed:  glob_prop_dict['feed']     = np.float(d_feed[0])
            if d_power: glob_prop_dict['power']    = np.float(d_power[0])
            if d_angle: glob_prop_dict['angle']    = np.float(d_angle[0])
            if d_split: glob_prop_dict['split']    = np.int(d_split[0])
            if d_radius:glob_prop_dict['radius'   ]= np.array([0,0,np.float(d_radius[0])])
            if d_cut_dir: glob_prop_dict['cut_dir']= d_cut_dir
            if d_coord: glob_prop_dict['ref_coord'] = [x for x in shape.insert]

            if d_key:
                forced_key = d_key[0]

    return glob_prop_dict, forced_key



def find_l_el(read_dir, el_kt_list, sorted_knots, master_knot):
    # find all elements including master_knot and put into el_index list
    el_index = [i for i, element in enumerate(
        el_kt_list) if master_knot in element]

    seg1, seg2 = elements_knots2coords(
        [el_kt_list[i] for i in el_index], sorted_knots)

    if cw_order(seg1, seg2) == read_dir:
        cur_ind = el_index[1]
    else:
        cur_ind = el_index[0]

    last_el = el_kt_list.pop(cur_ind)
    excl_knot = [x for x in last_el if x != master_knot]  # take the other knot

    return (last_el, excl_knot)


def cw_order(seg1, seg2):
    common_el = [x for x in list(set(seg1) & set(seg2))]
    u = sub_points(common_el, list(set(seg1) - set(common_el)))
    v = sub_points(common_el, list(set(seg2) - set(common_el)))
    if np.linalg.norm(np.cross(u, v)) > 0:
        return False
    else:
        return True

def merge(dict1, dict2):
    ''' Return a new dictionary by merging two dictionaries recursively. '''

    result = deepcopy(dict1)

    for key, value in dict2.items():
        if isinstance(value, collections.Mapping):
            result[key] = merge(result.get(key, {}), value)
        else:
            result[key] = deepcopy(dict2[key])

    return result
def update_dict(d, u):
    for k, v in u.items():
        if isinstance(v, collections.Mapping):
            default = v.copy()
            default.clear()
            r = update_dict(d.get(k, default), v)
            d[k] = r
        else:
            d[k] = v
    return d

def print_nested(val, buf, nesting = -1, prefix ='', bt=''):
    if len(val):
        nesting += 1
        for i, k in enumerate(sorted(val)):
            # print('val ',k)
            text_path= k
            if nesting==0:
                bt=' '*(nesting-1)
            else:
                if i+1 == len(val):
                    bt=prefix +'*--'
                else:
                    bt=prefix + '|--'
            buf.append('{}{}{}'.format(bt, text_path, ''))
            if len(val) > 1 and i != len(val) - 1:
                tmp_prefix = prefix + '| '
            else:
                tmp_prefix = prefix + ' '
            buf = print_nested(val[k], buf, nesting, tmp_prefix, bt)
    else:
        buf.append(val)
    return buf


def layers2seq(fname, req_layer):
    '''layer naming convention:
        A#XX#Y#ZZ~comment
        A - layername
        XX- sequence_number
        Y - column number 0 - flat, 1-XY, 2-UV
        ZZ- section number, might be empty
        ~comment
    '''

    key = r'(^{})#(\d+)#([01_])#((\d+)?)(.*)'.format(req_layer)

    layer_list = []
    seq_layer_dict = nested_dict()

    dwg = ezdxf.readfile(fname)
    for layer in dwg.layers:
        split_layer_name = re.findall(key, layer.dxf.name, re.IGNORECASE)
        if split_layer_name:
            seq_idx = split_layer_name[0][1]
            col_idx = split_layer_name[0][2]
            seq_layer_dict[seq_idx][col_idx][layer.dxf.name]=[]

    seq_layer_dict = OrderedDict(sorted(seq_layer_dict.items(), key = lambda x: x[0]))

    print('{:-^79}'.format('MODEL STRUCTURE'))
    for var in print_nested(seq_layer_dict, []):
        if var: print(var)

    seq_list = []
    for seq_key in sorted(seq_layer_dict.keys()):
        col_list = []
        for plane_key in sorted(seq_layer_dict[seq_key].keys()):
            layer_list = []
            for layer in sorted(seq_layer_dict[seq_key][plane_key].keys()):
                layer_list.append(layer)
            col_list.append(layer_list)
        seq_list.append(col_list)
    return seq_list, seq_layer_dict



if __name__ == '__main__':

    test_data = np.array([[[-1, -4, 0],[-1, -4, 1], [-2, -2, 2], [-1, -1, 3]],
                          [[4, -1, 0],[4, -1, 1], [2, -1, 2], [1, 0, 3]],
                          [[-2, 2, 0],[-2, 2, 1], [0, 3, 2], [0, 1, 3]]])

    buff_2ang, buff_2R, buff_2PS, buff_2v = transform(test_data, add_pedestal_bottom=True)
