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
  orig
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

    def origin(self, opt_list):

        stat_X= [min(self.data, key= lambda x: x[0])[0], max(self.data, key= lambda x: x[0])[0]]
        stat_Y= [min(self.data, key= lambda x: x[1])[1], max(self.data, key= lambda x: x[1])[1]]
        print stat_X
        print stat_Y
        offset_X = 0
        offset_Y = 0

        for opt in opt_list:
            if opt=='maxX': offset_X= -stat_X[1]
            if opt=='minX': offset_X= -stat_X[0]
            if opt=='maxY': offset_Y= -stat_Y[1]
            if opt=='minY': offset_Y= -stat_Y[0]
            if opt=='symX': offset_X= -(0.5 * sum(stat_X))
            if opt=='symY': offset_Y= -(0.5 * sum(stat_Y))

        self.translate([offset_X, offset_Y, 0])

#*****************************************parser = argparse.ArgumentParser(description='test')
parser = argparse.ArgumentParser(description='knots modifier')
parser.add_argument('-i', '--input', nargs='+',
                    type=str, help='input filename')
parser.add_argument('-o', '--output', nargs='+',
                    type=str, help='output filename. if not specified, the prog overrides the input file')
parser.add_argument('-s', '--scale', nargs='+', type=float, help='scale params: X Y Z')
parser.add_argument('-r', '--rotate', nargs='+', type=float, help='rotate params: X Y ang[deg]')
parser.add_argument('-t', '--translate', nargs='+',
                    type=float, help='translate params: X Y Z')
parser.add_argument('-org', '--origin', nargs='+',
                    type=str, help='adjust origin params (max 2 args): maxX|minX|maxY|minY|symX|symY')
args = parser.parse_args()

input_file_list = args.input
output_file = args.output
scale_args = args.scale
rotate_args = args.rotate
translate_args = args.translate
org_args = args.origin

if not input_file_list:
    dir_path = os.getcwd()
    input_file_list = [i for i in os.listdir(dir_path) if i.endswith('.knt')]

for input_file in input_file_list:
    dataxy = []
    print('{1}\nprocessing: {0}\n{1}'.format(input_file,'-'*24))
    dataxy = read_data(input_file)

    if dataxy:
        pathxy = path_points(dataxy)

        if scale_args and len(scale_args) == 3:
            pathxy.scale(scale_args)
            print('scale DONE')

        if rotate_args and len(rotate_args) == 3:
            pathxy.rotate(rotate_args)
            print('rotate DONE')

        if translate_args and len(translate_args) == 3:
            pathxy.translate(translate_args)
            print('translate DONE')

        if org_args:
            pathxy.origin(org_args)
            print('orig DONE')

        if not output_file:
            output_file = input_file

        write_data(output_file[0], pathxy.data, msg=True)

    output_file=[]
