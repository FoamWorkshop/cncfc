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


def knots2gcode(ct_path, name='gcode', global_header='False', subset_header='False'):
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
        for var in ct_path:
            f.write('G1 X{0:8.1f} Y{1:8.1f} U{0:8.1f} V{1:8.1f}\n'.format(
                var[0], var[1]))
    ##############################################
        if subset_header:
            f.write("\no<{0}> endsub\n".format(name))
            f.write("M2 (program end)")

        if global_header:
            f.write("M2 (program end)")

#*********************************************************************DEFAULT PARAMETERS
dflt_dxf_list = 'all'  # decimal accuracy
dflt_dec_acc = 3  # decimal accuracy
dflt_n_arc = 10  # number of segments
dflt_l_arc = 1  # minimal segment length
dflt_path_dir = 1  # closed path collecting direction
#*********************************************************************PROGRAM
knt_data = []
parser = argparse.ArgumentParser(description='test')
parser.add_argument('-i', '--input', nargs='+', help='input filenames')
parser.add_argument('-sh', '--subset_header', action='store_true')
parser.add_argument('-gh', '--global_header', action='store_true')

args = parser.parse_args()

knt_list = args.input
subset_header = args.subset_header
global_header = args.global_header

for knt_set in knt_list:
    print 'processing:', knt_set

    knt_data = read_data(knt_set, True)

    if knt_data:
        knots2gcode(knt_data,knt_set.replace('.knt',''), global_header, subset_header)
