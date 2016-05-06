#!/usr/bin/python
__author__ = 'FoamWorkshop'
'''the purpose of the program is to generate 3d model for FreeCAD from 2 given sections
modknot.py:
 -i input knots data file
 -o output knots data file

 -op:
  move [u,v,z]
  rotate [u,v,z] ang
  scale [sc_x, sc_y, sc_z]
 -rp: reverse path

info:
 the tool transforms selected knot path. the program can perform only one operation at time.
'''

import os
import argparse

import sys
import numpy as np
from cncfclib import *


class path_points:
    '''simple'''

    def __init__(self, data):
        self.data = data

    def rotate(self, rt):
        rc_x, rc_y, angl_deg = rt
        angl_rad = np.deg2rad(angl_deg)
        sn = np.sin(angl_rad)
        cs = np.cos(angl_rad)
        self.data = [(cs * (line[0] - rc_x) - sn * (line[1] - rc_y) + rc_x,
                      sn * (line[0] - rc_x) + cs * (line[1] - rc_y) + rc_y,
                      line[2]) for line in self.data]

    def translate(self, tr):
        self.data = [(line[0] + tr[0], line[1] + tr[1], line[2] + tr[2])
                     for line in self.data]

    def scale(self, sc):
        self.data = [(line[0] * sc[0], line[1] * sc[1], line[2] * sc[2])
                     for line in self.data]

#*********************************************************************DEFAULlist = 'all'  #decimal accuracy
dflt_output_file = 3  # decimal accuracydfl,number of segments
dflt_rotate = 1  # minimal segment lengthelp  1 #closed path collecting direction
#*****************************************parser = argparse.ArgumentParser(description='test')
parser = argparse.ArgumentParser(description='knots modifier')
parser.add_argument('-i', '--input', nargs='+',
                    type=str, help='input filename')
parser.add_argument('-o', '--output', nargs='+',
                    type=str, help='output filename')
parser.add_argument('-s', '--scale', nargs='+', type=float, help='scale')
parser.add_argument('-r', '--rotate', nargs='+', type=float, help='rotate')
parser.add_argument('-t', '--translate', nargs='+',
                    type=float, help='translate')

args = parser.parse_args()


input_file_list = args.input
output_file = args.output
scale_args = args.scale
rotate_args = args.rotate
translate_args = args.translate

for input_file in input_file_list:
    dataxy = []
    print('{1}\nprocessing: {0}\n{1}'.format(input_file,'-'*24))
    dataxy = read_data(input_file)

    if dataxy:
        pathxy = path_points(dataxy)

        if scale_args and len(scale_args) == 3:
            print('scale: [{0:.2f} {1:.2f} {2:.2f}]'.format(
                scale_args[0], scale_args[1], scale_args[2]))
            pathxy.scale(scale_args)

        if rotate_args and len(rotate_args) == 3:
            print('rotate: [{0:.2f} {1:.2f}] angle {2:.2f}'.format(
                rotate_args[0], rotate_args[1], rotate_args[2]))
            pathxy.rotate(rotate_args)

        if translate_args and len(translate_args) == 3:
            print('translate: [{0:.2f} {1:.2f} {2:.2f}]'.format(
                translate_args[0], translate_args[1], translate_args[2]))
            pathxy.translate(translate_args)

        if not output_file:
            output_file = input_file

        write_data(output_file, pathxy.data, msg=True)

    output_file=[]
