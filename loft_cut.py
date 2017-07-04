#!/usr/bin/python
from __future__ import division
FREECADPATH = '/usr/lib/freecad/lib/'  # path to your FreeCAD.so or FreeCAD.dll file

import sys
sys.path.append('/home/adam/Documents/00.projects/02.python/cncfc')
sys.path.append(FREECADPATH)
import cncfclib

import os
import argparse
import sys
import dxfgrabber
import numpy as np
from numpy import NaN, pi
from cncfclib import *
import collections
from scipy import interpolate
np.set_printoptions(precision=1)
import FreeCAD
import Part
from scipy.interpolate import CubicSpline

pt=collections.namedtuple('pt',['x', 'y'])

def angl_conv(x):
    out = x

    if x<0:
        out = 2*pi + x

    return out

vangl=np.vectorize(angl_conv)

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


def coords2file(name,coords_XU, coords_YV):
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
        print i
        P1 = np.vstack([P1_X[:,i], P1_Y[:,i]]).T
        P2 = np.vstack([P2_X[:,i], P2_Y[:,i]]).T
        O = np.vstack([O_X[:,i], O_Y[:,i]]).T
        C = cross_point(P1, P2, O)
        C_r[:,i]=radius(C,O)[:,0]
        C_a[:,i]=(angle_test( C_S, O, C) + C_a_ref)[:,0]* 180/pi

    # print C_a
        # C_a[:,i] = C_a[:,i]#+C_a_ref
    print C_r
    print C_a
    # #
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
    print C_Z
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

def dxf_read(files, layer_name, tol):

    line_count = 0
    circle_count = 0
    p_set= []
    c_set= []

    path_offset = (0,0,0)

    for shape in dxf.entities:
        if shape.layer == layer_name:

            if shape.dxftype == 'MTEXT':
                if shape.raw_text == 'coord_0':
                    # print('path offset: {}'.format(shape.insert))
                    path_offset = tuple(round(x, tol) for x in shape.insert)
                    print('path offset: {}'.format(path_offset))

    for shape in dxf.entities:
        if shape.layer == layer_name:

            if shape.dxftype == 'LINE':
                line_count += 1
                p_st = shape.start
                p_en = shape.end
                p_set.append( pt(round(p_st[0]-path_offset[0],tol), round(p_st[1]-path_offset[1],tol)))
                p_set.append( pt(round(p_en[0]-path_offset[0],tol), round(p_en[1]-path_offset[1],tol)))

            if shape.dxftype == 'CIRCLE':
                circle_count += 1
                pt_O_list = shape.center
                c_set.append(pt(round(pt_O_list[0]-path_offset[0],tol), round(pt_O_list[1]-path_offset[1],tol)))

    return (p_set, c_set, [line_count, circle_count])

def model_data_check(layer_P_list,layer_O_list, layer_c_set_list):
    for i, (pt_P_list, pt_O_list, c_set) in enumerate(zip(layer_P_list,layer_O_list, layer_c_set_list)):
        print layer_name
        if len(c_set) == 0:
            print('section: {}'.format(i))
            print('no circle defining the start point')
            return 0
        elif len(c_set) >1:
            print('section: {}'.format(i))
            print('more than one circle defining the start point:\n {}'.format(c_set))
            return 0
        elif not(c_set[0] in pt_P_list):
            print('section: {}'.format(i))
            print('circle center position does not match any section node')
            print('c_set: {}'.format(c_set[0]))
            print('P list:')
            for i, var in enumerate(pt_P_list):
                print('{}: {}'.format(i,var))

            return 0
        elif len(pt_O_list) != 1:
            print('section: {}'.format(i))
            print('error, not single center point: {}'.format(pt_O_list))
            return 0
        elif len(p_set)/2 != len(pt_P_list):
            print('section: {}'.format(i))
            print('some extra sections are present')
            return 0
        else:
            return 1
#*********************************************************************DEFAULT PARAMETERS
dflt_dxf_list = 'all'
dflt_layer_list= 'all'

dflt_dec_acc = 4  # decimal accuracy
dflt_n_arc = 10  # number of segments
dflt_l_arc = 0.1  # minimal segment length
dflt_path_dir = 1  # closed path collecting direction
dflt_sect=1
dflt_interp='linear'

prof_segment_list = []
prof_R_list = []
prof_th_list = []
start_offset_list = []
cross_point_list =[]
spoke_prof_list = []
C_list = []
C_a_list = []
C_r_list = []
O_list = []
layer_P_list=[]
layer_O_list=[]
layer_c_set_list=[]
seg_P1_list=[]
seg_P1_X_list=[]
seg_P1_Y_list=[]
seg_O_X_list=[]
seg_O_Y_list=[]
#*********************************************************************PROGRAM
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                description='''\
--------------------------
lofted shape rules for dxf
--------------------------

the lofted shape is defined by:

1. profile layers:
[prefix_XX_YYYY] where: XX   - part numeber; YYYY - profile number.
Each profile is build with n LINES and a CIRCLE. The other shapes are neglected. he profile is builf with n LINES. The lines have one end coincident (O_POINT) and the other is coincident with the cut contour of the profile. There is one distiguinshed line (0_LINE) which end is marked with a CIRCLE center. This line is treated as the first and the other lines are ordered counterclockwise. All profiles must be build with THE SAME NUMBER of LINES.

      2. spin layer:
[prefix_XX_spin] where: XX   - part numeber
The spin is built with n CIRCLES which define JOINTS. The other shapes are neglected. Y position of the CIRCLE center points define coresponding positions of the profile O_POINTS. O_POINTs are assigned to JOINTS as the ordered profile names (LAYER NAMES) along the Y axis. Number of joints must correspond to the number of the profiles.

''',
                                usage='%(prog)s [-i dxf_filename] [-l lofted_shape_prefix]')



parser.add_argument('-i', '--input', type=str, required=True, help='input filename')
parser.add_argument('-l', '--profile_prefix', type=str, required=True, help='profile data prefix. prefix_xxxx - profile; prefix_spin - loft spin data')
parser.add_argument('-a', '--accuracy', type=int, default=dflt_dec_acc, help='decimal accuracy, default: 3')
parser.add_argument('-sect', '--sect', type=int, default=dflt_sect, help='number of additional interpolation sections along the spin, default: 0')
parser.add_argument('-interp', '--interp', type=str, default=dflt_interp, help='cutting path interpolation method between profiles: linear (default), nearest, zero, slinear, quadratic, cubic')

args = parser.parse_args()

dxf_list = args.input
layer_list = args.profile_prefix
dec_acc = args.accuracy
n_sect=args.sect
interp_meth=args.interp

dir_path = os.getcwd()
dxf_files = [i for i in os.listdir(dir_path) if i.endswith('.dxf')]

print(dxf_list)
if not dxf_list:
    print('the input file is not specified, use: -i "file name"')
# elif not (dxf_list in dxf_files):
#     print('the current dir does not include requested dxf files')
# elif not layer_list:
#     print('the profile prefix is not specified, use: -l "profile prefix"')



else:
    prof_pref= layer_list
    prof_spin= layer_list+'spin'

    files_dxf_member = dxf_list

    print("\n{0:24} {1:8}".format('extracting from: ',files_dxf_member))
    print('{0}'.format('-' * 13*5))

    case_name = os.path.splitext(files_dxf_member)
    dxf = dxfgrabber.readfile(files_dxf_member, {"assure_3d_coords": True})
    dxf_layers = dxf.layers

    # prof_layer_name_list = [var.name for var in dxf_layers if prof_pref in var.name ]

    prof_layer_name_list = [var.name for var in dxf_layers if prof_pref in var.name and var.name.replace(prof_pref,'').isnumeric()]
    prof_spin_layer_name_list = [var.name for var in dxf_layers if prof_spin == var.name]

    if not prof_layer_name_list:
        print('dxf file does not include specified profile layers: {}*'.format(prof_pref))

    elif len(prof_layer_name_list)<2:
        print('alt least 2 profiles must be defined')

    elif not prof_spin_layer_name_list:
        print('dxf file does not include specified spin layer: {}_spin'.format(prof_spin))

    elif len(prof_spin_layer_name_list)>1:
        print('dxf file includes more than 1 spin layer: {}'.format(prof_spin_layer_name_list))

    else:
        print('profile layers: {}'.format(sorted(prof_layer_name_list)))
        print('spin_layer: {}'.format(prof_spin_layer_name_list))

        for layer_name in sorted(prof_layer_name_list):
            p_set, c_set, shape_count = dxf_read(dxf, layer_name, dec_acc)
            #find ends of spokes P[0] - P[n]
            pt_P_list = [x for x, y in collections.Counter(p_set).items() if y == 1]
            #find center point 0
            pt_O_list = [x for x, y in collections.Counter(p_set).items() if y >  1]

            layer_P_list.append(pt_P_list)
            layer_O_list.append(pt_O_list)
            layer_c_set_list.append(c_set)

        #-----SEKCJA DO POPRAWY. WCZYTANIE KOORDYNATOW Z
        for layer_name in sorted(prof_spin_layer_name_list):
            dummy_1, c_set, dummy_2 = dxf_read(dxf, layer_name, dec_acc)
            C_Z = np.sort([var.y for var in c_set])
        #-----SEKCJA DO POPRAWY. WCZYTANIE KOORDYNATOW Z

        if model_data_check(layer_P_list,layer_O_list, layer_c_set_list):

            for pt_P_list, pt_O_list, c_set in zip(layer_P_list,layer_O_list, layer_c_set_list):

                n_spokes = len(pt_P_list)
                O=np.ones((n_spokes , 2)) * np.array(pt_O_list)
                S=np.ones((n_spokes , 2)) * np.array(c_set)
                P=np.ones((n_spokes , 2)) * np.array(pt_P_list)

                r = radius(P, O)
                a = angle_test( S, O, P)*180/pi
                P_ref = np.array([[1,0]])
                a_ref = angle_test( P_ref, O, S) * 180/pi
                a_inds = np.argsort(a.T)
                seg_P1 = P[a_inds][0]
                seg_P2 = np.roll(seg_P1,-1, axis=0)

                seg_P1_X_list.append(seg_P1[:,0])
                seg_P1_Y_list.append(seg_P1[:,1])

                seg_O_X_list.append(np.array(pt_O_list[0])[0])
                seg_O_Y_list.append(np.array(pt_O_list[0])[1])

                C = cross_point(seg_P1, seg_P2, O)
                C_S = np.ones((n_spokes , 2)) * C[0,:]
                C_a = angle_test( C_S, O, C) * 180/pi
                C_a_ref = angle_test( P_ref, O, C_S) * 180/pi
                C_r = radius(C,O)


                C_r_list.append(C_r.T)
                C_a_list.append(C_a.T + C_a_ref.T)

            # print seg_P1_list
            # print(np.vstack(seg_P1_X_list))
            # print(np.vstack(seg_P1_Y_list))

            interp_points(seg_P1_X_list, seg_P1_Y_list , seg_O_X_list, seg_O_Y_list, C_Z, n_sect, interp_meth,files_dxf_member,layer_list)

            interp_points_fc(seg_P1_X_list, seg_P1_Y_list, C_Z, n_sect, interp_meth,files_dxf_member)
            # C_a1=np.vstack(C_a_list)
            # # C_a2=np.roll(C_a1,-1, axis=0)
            #
            # C_r1=np.vstack(C_r_list)
            # # C_r2=np.roll(C_r1,-1, axis=0)
            #
            # yv = np.hstack((np.linspace(C_Z[0],C_Z[-1],n_sect),C_Z))
            # yv = np.unique(yv)
            # yv = np.sort(yv)
            #
            # print n_spokes
            # cut_in_swing = True
            # for i in range(n_spokes):
            #     print len(C_Z)
            #     print len(C_r1[:,i])
            #     path_C_r1 = interpolate.interp1d(C_Z, C_r1[:,i],kind=interp_meth)(yv)
            #     path_C_a1 = interpolate.interp1d(C_Z, C_a1[:,i],kind=interp_meth)(yv)
            #     f_name_C_r1='{0}_xyuv_{1:{fill}{align}4}.knt'.format(layer_list,i,fill='0',align='>')
            #     f_name_C_a1='{0}_b_{1:{fill}{align}4}.knt'.format(layer_list,i,fill='0',align='>')
            #     print('yv: {}'.format(yv))
            #     print('r:  {}'.format(path_C_r1))
            #     print('a:  {}'.format(path_C_a1))
            #
            #     if cut_in_swing and i%2:
            #         coords2file(f_name_C_r1, np.flipud(path_C_r1), np.flipud(yv))
            #         Rcoords2file(f_name_C_a1, np.flipud(path_C_a1))
            #     else:
            #         coords2file(f_name_C_r1, path_C_r1, yv)
            #         Rcoords2file(f_name_C_a1, path_C_a1)

print "\n end of program. thank you!"
