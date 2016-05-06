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
            f.write('G1 X{0:8.1f} Y{1:8.1f} U{2:8.1f} V{3:8.1f}\n'.format(
                varxy[0], varxy[1], varuv[0], varuv[1]))
    ##############################################
        if subset_header:
            f.write("\no<{0}> endsub\n".format(name))
            f.write("M2 (program end)")

        if global_header:
            f.write("M2 (program end)")




def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = x_limits[1] - x_limits[0]; x_mean = np.mean(x_limits)
    y_range = y_limits[1] - y_limits[0]; y_mean = np.mean(y_limits)
    z_range = z_limits[1] - z_limits[0]; z_mean = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinityl
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_mean - plot_radius, x_mean + plot_radius])
    ax.set_ylim3d([y_mean - plot_radius, y_mean + plot_radius])
    ax.set_zlim3d([z_mean - plot_radius, z_mean + plot_radius])

def plot_path3d(data1, data2, data3, data4):
    x1=[x for [x, y, z] in data1]
    y1=[y for [x, y, z] in data1]
    z1=[z for [x, y, z] in data1]

    x2=[x for [x, y, z] in data2]
    y2=[y for [x, y, z] in data2]
    z2=[z for [x, y, z] in data2]

    mx1=[x for [x, y, z] in data3]
    my1=[y for [x, y, z] in data3]
    mz1=[z for [x, y, z] in data3]

    mx2=[x for [x, y, z] in data4]
    my2=[y for [x, y, z] in data4]
    mz2=[z for [x, y, z] in data4]

    fig=plt.figure()
    ax=fig.gca(projection='3d')
    pltxy=ax.plot(x1, y1, z1,label='xy')
    pltuv=ax.plot(x2, y2, z2,label='uv')
    mashpltxy=ax.plot(mx1, my1, mz1,label='mxy')
    mashpltuv=ax.plot(mx2, my2, mz2,label='muv')

    ax.legend()
    set_axes_equal(ax)
    plt.show()
    plt.hold(True)

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
            print l0, l1
            tmp.append(p_l_intersection(p0,vec_n,l0,l1))
        return tmp
    else:
        return [0,0,0]

def gcodexyuv(dataxy, datauv):

    if len(dataxy)==len(datauv):
        print "data ok"
        tmp=[]
        for i in range(len(dataxy)):
            tmp.append('g1 x{0:6.2f} y{1:6.2f} u{2:6.2f} v{3:6.2f}'.format(dataxy[i][0], dataxy[i][1], datauv[i][0], datauv[i][1]))

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
d=500
d_rat = 0.5
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

args = parser.parse_args()

knt_list = args.input

subset_header = args.subset_header
global_header = args.global_header
center_model = args.center_model
output_f_name = args.output

d = args.distance
d_rat = args.distance_ratio


if len(knt_list)==1:
    knt_set_xy = knt_list[0]
    knt_set_uv = knt_list[0]

elif len(knt_list)>=2:
    knt_set_xy = knt_list[0]
    knt_set_uv = knt_list[1]

if not output_f_name:
    output_f_name='_'.join([knt_set_xy.replace('.knt',''),knt_set_uv.replace('.knt','')])

knt_data_xy = read_data(knt_set_xy, True)
knt_data_uv = read_data(knt_set_uv, True)

if len(knt_set_xy)!=len(knt_set_uv):
    print('knots: {0} - {1} are not balanced. EXIT'.format(knt_set_xy, knt_set_uv))
else:
    print('processing knots: {0} - {1}'.format(knt_set_xy, knt_set_uv))

    if center_model:
        # for varxy, varuv in zip(knt_data_xy, knt_data_uv):
        pool=zip(knt_data_xy, knt_data_uv)
        knt_data_xy=[[varxy[0], varxy[1], varxy[2]-0.5*(varxy[2]+varuv[2])] for varxy, varuv in pool]
        knt_data_uv=[[varuv[0], varuv[1], varuv[2]-0.5*(varxy[2]+varuv[2])] for varxy, varuv in pool]


    mashpathxy=p_l_intersection_series([0,0,d *  d_rat   ],[0,0,1],knt_data_xy,knt_data_uv)
    mashpathuv=p_l_intersection_series([0,0,d * (d_rat-1)],[0,0,1],knt_data_uv,knt_data_xy)

    write_data(knt_set_xy.replace('.knt','_xy.knt'),mashpathxy,True)
    write_data(knt_set_uv.replace('.knt','_uv.knt'),mashpathuv,True)
    knots2gcode(mashpathxy, mashpathuv, output_f_name, global_header, subset_header)
