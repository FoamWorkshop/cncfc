#!/usr/bin/python
__author__ = 'FoamWorkshop'

'''dxf2path.py, 6.04.2016 author: Adam Narozniak
dxf2path ptogram is a part of the CNCFCgcode generator. It automaticaly extracts cutting path from a dxf file.
The cutting path is split into:
1. IO_path - in/out path begining with single knot
2. ct_path - closed loop, begin and end in master knot position
the output of the program is a set of files with an ordered list of knots'''

import dxfgrabber
import numpy as np
import argparse
from cncfclib import *


def cross_prod(u, v):
    return u[0] * v[1] - u[1] * v[0]


def sub_points(p1, p2):
    vect = []
    p1 = [x for x in p1[0]]
    p2 = [x for x in p2[0]]
 #   print p3, p4
    print len(p1)
    print len(p2)
    if len(p1) == len(p2):
        for i, n in enumerate(p2):
            vect.append(n - p1[i])
        return vect
    return len(p1) * [None]


def knots2gcode(ct_pathxy, ct_pathuv, name='gcode', global_header='False', subset_header='False'):
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
        for varxy, varuv in zip(ct_pathxy, ct_pathuv):
            f.write('G1 X{0:8.3f} Y{1:8.3f} U{2:8.3f} V{3:8.3f}\n'.format(
                varxy[0], varxy[1], varuv[0], varuv[1]))
    ##############################################
        if subset_header:
            f.write("\no<{0}> endsub\n".format(name))
            f.write("M2 (program end)")

        if global_header:
            f.write("M2 (program end)")



def p_l_intersection(p0,vec_n,l0,l1):
    vec_l=np.subtract(l1,l0)
    param1=np.subtract(p0,l0)
    d=(np.dot(param1,vec_n))/(np.dot(vec_l,vec_n))
    vec_l=np.multiply(d,vec_l)
    return np.add(vec_l,l0)

def p_l_intersection_series(p0,vec_n,data1,data2):
    if len(data1)==len(data2):
        print "data ok"
        tmp=[]
        for i in range(len(data1)):
            l0=data1[i]
            l1=data2[i]
        #    print l0, l1
            tmp.append(p_l_intersection(p0,vec_n,l0,l1))
        return tmp
    else:
        return [0,0,0]

def gcodexyuv(dataxy, datauv):

    if len(dataxy)==len(datauv):
        print "data ok"
        tmp=[]
        for i in range(len(dataxy)):
            tmp.append('g1 x{0:6.3f} y{1:6.3f} u{2:6.3f} v{3:6.3f}'.format(dataxy[i][0], dataxy[i][1], datauv[i][0], datauv[i][1]))

        fgcode=open('test.ngc','w')
        for line in tmp:
            fgcode.write(line)
            fgcode.write('\n')
        fgcode.close()
        return tmp

    else:
        print "nie mozna wygenerowac g codu. rozne dlugosci sciezek."
        return 0

#*********************************************************************DEFAULT PARAMETERS
dflt_dxf_list = 'all'  # decimal accuracy
dflt_dec_acc = 3  # decimal accuracy
dflt_n_arc = 10  # number of segments
dflt_l_arc = 1  # minimal segment length
dflt_path_dir = 1  # closed path collecting direction
d=421
d_rat = 0.5
symm_pref = 's'
#*********************************************************************PROGRAM
knt_data = []
parser = argparse.ArgumentParser(description='test')
parser.add_argument('-i', '--input', nargs='+',type=str, help='input filenames')
parser.add_argument('-o', '--output', type=str, help='input filenames')
parser.add_argument('-sh', '--subset_header', action='store_true')
parser.add_argument('-gh', '--global_header', action='store_true')
parser.add_argument('-cm', '--center_model', action='store_true')
parser.add_argument('-d', '--distance', type=float, default=d, help='distance between columns')
parser.add_argument('-dr', '--distance_ratio', type=float, default=d_rat, help='(xy-C)/d')
parser.add_argument('-symm', '--symmetry', action='store_true')

args = parser.parse_args()

knt_list = args.input

subset_header = args.subset_header
global_header = args.global_header
center_model = args.center_model
output_f_name = args.output
symm_stat = args.symmetry

d = args.distance
d_rat = args.distance_ratio


if len(knt_list)==1:
    knt_set_xy = knt_list[0]
    knt_set_uv = knt_list[0]

elif len(knt_list)>=2:
    knt_set_xy = knt_list[0]
    knt_set_uv = knt_list[1]
    print '1xy:',knt_set_xy
    print '1uv:',knt_set_uv

if not output_f_name:
    output_f_name='_'.join([knt_set_xy.replace('.knt',''),knt_set_uv.replace('.knt','')])

knt_data_xy = read_data(knt_set_xy, False)
knt_data_uv = read_data(knt_set_uv, False)

if len(knt_set_xy)!=len(knt_set_uv):
    print('knots: {0} - {1} are not balanced. EXIT'.format(knt_set_xy, knt_set_uv))
else:
    print('processing knots: {0} - {1}'.format(knt_set_xy, knt_set_uv))

    pool=zip(knt_data_xy, knt_data_uv)
    knt_data_xy=[[varxy[0], varxy[1], varxy[2]+int(not(varxy[2]-varuv[2]))] for varxy, varuv in pool]
    knt_data_uv=[[varuv[0], varuv[1], varuv[2]-int(not(varxy[2]-varuv[2]))] for varxy, varuv in pool]

    pool=[]
    if center_model:
        # for varxy, varuv in zip(knt_data_xy, knt_data_uv):
        pool=zip(knt_data_xy, knt_data_uv)
        knt_data_xy=[[varxy[0], varxy[1], varxy[2]-0.5*(varxy[2]+varuv[2])] for varxy, varuv in pool]
        knt_data_uv=[[varuv[0], varuv[1], varuv[2]-0.5*(varxy[2]+varuv[2])] for varxy, varuv in pool]

    mashpathxy=p_l_intersection_series([0,0,d *  d_rat   ],[0,0,1],knt_data_xy,knt_data_uv)
    mashpathuv=p_l_intersection_series([0,0,d * (d_rat-1)],[0,0,1],knt_data_uv,knt_data_xy)

    write_data('{0}_{1}.knt'.format('xy',output_f_name),mashpathxy,True)
    write_data('{0}_{1}.knt'.format('uv',output_f_name),mashpathuv,True)

    knots2gcode(mashpathxy, mashpathuv, output_f_name, global_header, subset_header)

    if symm_stat:
        knots2gcode(mashpathuv, mashpathxy, ''.join([symm_pref, output_f_name]), global_header, subset_header)
