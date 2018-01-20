import sys
import os
import numpy as np
from stl import mesh
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from collections import Iterable

def flatten(items):
    """Yield items from any nested iterable; see REF."""
    for x in items:
        if isinstance(x, Iterable) and not isinstance(x, (str, bytes)):
            yield from flatten(x)
        else:
            yield x

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
    '''vector to plane projection.\n
    input: u - vector, n - plane normal, O - attachement point\n
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

def add_pedestal(pos, add_pedestal_top=True, add_pedestal_bottom=True):
    r = 35
    h = np.array([0, 5, 5])
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




if __name__ == '__main__':

    # cross_point(P1, P2, Q1, Q2)
    test_data = np.array([[[-1, -4, 0],[-1, -4, 1], [-2, -2, 2], [-1, -1, 3]],
                          [[4, -1, 0],[4, -1, 1], [2, -1, 2], [1, 0, 3]],
                          [[-2, 2, 0],[-2, 2, 1], [0, 3, 2], [0, 1, 3]]])
    # test_data = np.array([[[-2, -2, 2], [-1, -1, 3]],
    #                       [[2, -1, 2], [1, 0, 3]],
    #                       [[0, 3, 2], [0, 1, 3]]])
    # print(test_data)

    buff_2ang, buff_2R, buff_2PS, buff_2v = transform(test_data, add_pedestal_bottom=True)

    # print('ang', buff_2ang)
    # print('buf2', buff_2R)
    # print('z', buff_2PS)
    # print('2v', buff_2v  )
