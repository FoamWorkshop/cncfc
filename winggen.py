#!/usr/bin/python
from __future__ import division
'''the purpose of the program is to generate 3d model for FreeCAD from 2 given sections'''
FREECADPATH = '/usr/lib/freecad/lib/'  # path to your FreeCAD.so or FreeCAD.dll file
import sys
sys.path.append(FREECADPATH)
import os
import numpy as np
import numpy.linalg as npl
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt

import FreeCAD
import Part
from cncfclib import *
import argparse
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import interp1d
from scipy.spatial import distance


def set_axes_equal(ax):
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = x_limits[1] - x_limits[0]
    x_mean = np.mean(x_limits)
    y_range = y_limits[1] - y_limits[0]
    y_mean = np.mean(y_limits)
    z_range = z_limits[1] - z_limits[0]
    z_mean = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinityl
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_mean - plot_radius, x_mean + plot_radius])
    ax.set_ylim3d([y_mean - plot_radius, y_mean + plot_radius])
    ax.set_zlim3d([z_mean - plot_radius, z_mean + plot_radius])


def plot_path2d(data1, data2, data3, data4):
    x1 = [x for [x, y, z] in data1]
    y1 = [y for [x, y, z] in data1]
    z1 = [z for [x, y, z] in data1]

    x2 = [x for [x, y, z] in data2]
    y2 = [y for [x, y, z] in data2]
    z2 = [z for [x, y, z] in data2]

    mx1 = [x for [x, y, z] in data3]
    my1 = [y for [x, y, z] in data3]
    mz1 = [z for [x, y, z] in data3]

    mx2 = [x for [x, y, z] in data4]
    my2 = [y for [x, y, z] in data4]
    mz2 = [z for [x, y, z] in data4]

    pltxy = plt.plot(x1, y1, 'ro-', label='xy')
    pltuv = plt.plot(x2, y2, 'bs-', label='uv')
    mpltxy = plt.plot(mx1, my1, 'go-', label='cut_xy')
    mpltuv = plt.plot(mx2, my2, 'ks-', label='cut_uv')

    plt.legend()
    plt.axis('equal')
    plt.axis([min([min(x1), min(x2), min(mx1), min(mx2)]),
              max([max(x1), max(x2), max(mx1), max(mx2)]),
              min([min(y1), min(y2), min(my1), min(my2)]),
              max([max(y1), max(y2), max(my1), max(my2)])])
    plt.grid(True)
    plt.interactive(True)
    plt.show(block=False)
    plt.hold(True)


def plot_path3d(data1, data2, data3, data4):
    x1 = [x for [x, y, z] in data1]
    y1 = [y for [x, y, z] in data1]
    z1 = [z for [x, y, z] in data1]

    x2 = [x for [x, y, z] in data2]
    y2 = [y for [x, y, z] in data2]
    z2 = [z for [x, y, z] in data2]

    mx1 = [x for [x, y, z] in data3]
    my1 = [y for [x, y, z] in data3]
    mz1 = [z for [x, y, z] in data3]

    mx2 = [x for [x, y, z] in data4]
    my2 = [y for [x, y, z] in data4]
    mz2 = [z for [x, y, z] in data4]

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    pltxy = ax.plot(x1, y1, z1, 'ro-', label='xy')
    pltuv = ax.plot(x2, y2, z2, 'bs-', label='uv')
    mashpltxy = ax.plot(mx1, my1, mz1, label='mxy')
    mashpltuv = ax.plot(mx2, my2, mz2, label='muv')

    ax.legend()
    set_axes_equal(ax)
    plt.show()
    plt.hold(True)


def p_l_intersection(p0, vec_n, l0, l1):
    vec_l = np.subtract(l1, l0)
    param1 = np.subtract(p0, l0)
    d = (np.dot(param1, vec_n)) / (np.dot(vec_l, vec_n))
    vec_l = np.multiply(d, vec_l)
    return np.add(vec_l, l0)


def p_l_intersection_series(p0, vec_n, data1, data2):
    if len(data1) == len(data2):
        print "data ok"
        tmp = []
        for i in range(len(data1)):
            l0 = data1[i]
            l1 = data2[i]
            print l0, l1
            tmp.append(p_l_intersection(p0, vec_n, l0, l1))
        return tmp
    else:
        return [0, 0, 0]


class model:

    def __init__(self, sect1, sect2):
        self.sect1 = sect1
        self.sect2 = sect2


class path_points:
    '''simple'''

    def __init__(self, data):
        self.data = data

    def rotate(self, rc_x, rc_y, angl_deg):
        angl_rad = np.deg2rad(angl_deg)
        sn = np.sin(angl_rad)
        cs = np.cos(angl_rad)
        self.data = [(cs * (line[0] - rc_x) - sn * (line[1] - rc_y) + rc_x,
                      sn * (line[0] - rc_x) + cs * (line[1] - rc_y) + rc_y,
                      line[2]) for line in self.data]

    def translate(self, tr_x, tr_y, tr_z):
        self.data = [(line[0] + tr_x, line[1] + tr_y, line[2] + tr_z)
                     for line in self.data]

    def scale(self, sc_x, sc_y):
        self.data = [(line[0] * sc_x, line[1] * sc_y, line[2])
                     for line in self.data]


def gcodexyuv(dataxy, datauv):

    if len(dataxy) == len(datauv):
        print "data ok"
        tmp = []
        for i in range(len(dataxy)):
            tmp.append('g1 x{0:6.2f} y{1:6.2f} u{2:6.2f} v{3:6.2f}'.format(
                dataxy[i][0], dataxy[i][1], datauv[i][0], datauv[i][1]))

        fgcode = open('test.ngc', 'w')
        for line in tmp:
            fgcode.write(line)
            fgcode.write('\n')
        fgcode.close()
        return tmp

    else:
        print "nie mozna wygenerowac g codu. rozne dlugosci sciezek."
        return 0


def PolyArea(x, y):
    return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


def silentremove(filename):
    if os.path.exists(filename):
        os.remove(filename)


def points2face(points):
    fc_profile1_elems = []

    for p1, p2 in zip(points, points[1:]):
        fc_profile1_elems.append(Part.makeLine(p1, p2))

    fc_profile1_wire = Part.Wire(fc_profile1_elems)
    return Part.Face(fc_profile1_wire)


def balance_points(data1, data2):
    '''function balances 2 data sets, data2 to data1'''
    dist_1 = []
    dist_2 = []
    buff = []
    x = []
    y = []
    z = []
    data1_u = []
    data1_N = len(data1)
    data2_N = len(data2)
    cos_uv = []
    #  length_set=[]
    dist_1.append(0)
    dist_1.extend([distance.euclidean(u, v) for u, v in zip(data1, data1[1:])])
    length_1 = np.cumsum(dist_1)

    dist_2.append(0)
    dist_2.extend([distance.euclidean(u, v) for u, v in zip(data2, data2[1:])])
    length_2 = np.cumsum(dist_2)
    # print data1[6][0]
    data1_A = 0

    for i in range(data1_N - 1):
        data1_A += 0.5 * (data1[i][0] * data1[i + 1]
                          [1] - data1[i][1] * data1[i + 1][0])

    data1_C1 = 0
    data1_C2 = 0
    for i in range(data1_N - 1):
        data1_C1 += 1 / (6 * data1_A) * (data1[i][0] + data1[i + 1][0]) * (
            data1[i][0] * data1[i + 1][1] - data1[i][1] * data1[i + 1][0])
        data1_C2 += 1 / (6 * data1_A) * (data1[i][1] + data1[i + 1][1]) * (
            data1[i][0] * data1[i + 1][1] - data1[i][1] * data1[i + 1][0])

    print("total length data1: {0},\narea: {1},\nC1: {2}\nC2: {3}".format(
        length_1[-1], data1_A, data1_C1, data1_C2))

    data1_u = [(var[0] - data1_C1, var[1] - data1_C2) for var in data1]
    data2_u = [(var[0] - data1_C1, var[1] - data1_C2) for var in data2]

    print("\n ====data1 vectors===")

    for var in data1_u:
        print var

    print("\n ====data2 vectors===")

    for var in data2_u:
        print var

    # data1_x=[ var[0] for var in data1]
    # data1_y=[ var[1] for var in data1]
    # data1_z=[ var[2] for var in data1]

   # print data1_x

parser = argparse.ArgumentParser(description='half wing generator')

parser.add_argument('-i', '--input_airfoils', type=str, nargs='+',
                    required=True, help='input root and tip airfoil knots')
parser.add_argument('-o', '--output_airfoils', type=str,
                    nargs='+', help='output root and tip airfoil knots')
parser.add_argument('-fco', '--fc_model', type=str, help='FreeCAD output model')
parser.add_argument('-cr', '--root_coord', type=float,
                    default=1, help='rood coord')
parser.add_argument('-ct', '--tip_coord', type=float,
                    default=1, help='tip coord')
parser.add_argument('-b', '--span', type=float, default=1, help='wing span')
parser.add_argument('-s', '--area', type=float, default=1, help='wing area')
parser.add_argument('-ar', '--aspect_ratio', type=float,
                    default=1, help='aspect ratio')
parser.add_argument('-lambda', '--taper_ratio', type=float,
                    default=1, help='aspect ratio')
parser.add_argument('-theta', '--twist_angle', type=float,
                    default=0, help='twist angle')
parser.add_argument('-sa', '--sweep_angle', type=float,
                    default=0, help='sweep angle')
parser.add_argument('-d', '--dihedral', type=float,
                    default=0, help='wing dihedral')
parser.add_argument('-dx', '--tip_dx', type=float,
                    default=0, help='sweep angle')
parser.add_argument('-dy', '--tip_dy', type=float,
                    default=0, help='wing dihedral')

args = parser.parse_args()

input_airfoils = args.input_airfoils
output_airfoils = args.output_airfoils
fc_model =args.fc_model

if len(input_airfoils) == 1:
    i_airfoil_root, i_airfoil_tip = input_airfoils[0], input_airfoils[0]
else:
    i_airfoil_root, i_airfoil_tip = input_airfoils[:2]

if output_airfoils:
    if len(output_airfoils) == 1:
        o_airfoil_root, o_airfoil_tip = output_airfoils[0], output_airfoils[0]
    else:
        o_airfoil_root, o_airfoil_tip = output_airfoils[:2]
else:
    o_airfoil_root = '_root'.join(i_airfoil_root.split('.'))
    o_airfoil_tip = '_tip'.join(i_airfoil_tip.split('.'))

if not fc_model:
    fc_model= '_'.join([i_airfoil_root.split('.')[0], i_airfoil_tip.split('.')[0]])

c_r = args.root_coord
c_t = args.tip_coord
b = args.span
s = args.area
ar = args.aspect_ratio
lambd = args.taper_ratio
theta = args.twist_angle
se = args.sweep_angle
gamma = args.dihedral
tip_dx = args.tip_dx
tip_dy = args.tip_dy

data_tip = read_data(i_airfoil_tip)
data_root = read_data(i_airfoil_root)

if len(data_tip) != len(data_root):
    print('number of root and tip airfoils do not match')
else:
    if data_tip[0]!=data_tip[-1]:
        print('tip profile is open')

    elif data_root[0]!=data_root[-1]:
        print('root profile is open')

    else:
        root_profile = path_points(data_root)
        tip_profile = path_points(data_tip)

        root_profile.scale(c_r, c_r)
        # root_profile.rotate(0, 0, 0)
        # root_profile.translate(0, 0, 0)

        tip_profile.scale(c_t, c_t)
        tip_profile.rotate(0, 0, theta)
        tip_profile.translate(tip_dx, tip_dy, b / 2)

        doc = FreeCAD.newDocument()

        myPart = doc.addObject('Part::Feature', 'wing')

        points1 = root_profile.data
        points2 = tip_profile.data
        print 'root points'
        for var in points1:
            print var
        print 'tip points'
        for var in points2:
            print var


        face1 = points2face(points1)
        face2 = points2face(points2)

        poly1 = Part.makePolygon(points1)
        poly2 = Part.makePolygon(points2)

        myPart.Shape = Part.makeLoft([poly1, poly2], True)
        # silentremove('test_data1')
        doc.saveAs(fc_model + '.fcstd')
