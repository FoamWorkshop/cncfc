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
import cncfclib as cf
import collections
from scipy import interpolate
np.set_printoptions(precision=1)
import FreeCAD
import Part
from scipy.interpolate import CubicSpline

pt=collections.namedtuple('pt',['x', 'y'])


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

                r = cf.radius(P, O)
                a = cf.angle_test( S, O, P)*180/pi
                P_ref = np.array([[1,0]])
                a_ref = cf.angle_test( P_ref, O, S) * 180/pi
                a_inds = np.argsort(a.T)
                seg_P1 = P[a_inds][0]
                seg_P2 = np.roll(seg_P1,-1, axis=0)

                seg_P1_X_list.append(seg_P1[:,0])
                seg_P1_Y_list.append(seg_P1[:,1])

                seg_O_X_list.append(np.array(pt_O_list[0])[0])
                seg_O_Y_list.append(np.array(pt_O_list[0])[1])

                C = cf.cross_point(seg_P1, seg_P2, O)
                C_S = np.ones((n_spokes , 2)) * C[0,:]
                C_a = cf.angle_test( C_S, O, C) * 180/pi
                C_a_ref = cf.angle_test( P_ref, O, C_S) * 180/pi
                C_r = radius(C,O)


                C_r_list.append(C_r.T)
                C_a_list.append(C_a.T + C_a_ref.T)

            cf.interp_points(seg_P1_X_list, seg_P1_Y_list , seg_O_X_list, seg_O_Y_list, C_Z, n_sect, interp_meth,files_dxf_member,layer_list)
            cf.interp_points_fc(seg_P1_X_list, seg_P1_Y_list, C_Z, n_sect, interp_meth,files_dxf_member)

print "\n end of program. thank you!"
