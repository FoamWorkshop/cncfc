#!/usr/bin/python

from __future__ import division
__author__ = 'FoamWorkshop'
'''
program options:
should be added to manipulation script?

-i [input file[ input file]] - top view profile, optional side view profile.
                               if the side view profile is a wire and z=0, the surf is treated as a simple fold
                               if the side view profile is as wire and z!=0, the surf is treated as a ruled surface

-o [output file] - unfolded shape

freecad options:
    -fc [part name] - create freecad output with specified part name, if not specified, the part name is top profile
    -fct [sheet thickness] - default 0 creates plane, thickness side as 'the right hand rule' of the path algorithm

TODO:
    -add an output with a crise pattern
    -merge with 'knots manipulation script'
    -rename variables, make comments
    -organize functionalities in a user friendly way
    -check if required modules are installed and disable/enable appropriate functionalities
'''

# path to your FreeCAD.so or FreeCAD.dll file
FREECADPATH = '/usr/lib/freecad/lib/'
import sys
sys.path.append(FREECADPATH)
import FreeCAD
from FreeCAD import Base
import Part
import os
from cncfclib import *
import argparse
import dxfgrabber
import numpy as np


def cross_prod(u, v):
    return u[0] * v[1] - u[1] * v[0]


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
            knots_rank.append(
                [var[0], [x for a in el_kt_list for x in a].count(var[0])])
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
        f.write('{0:.2f} {1:.2f}\n'.format(coord[0], coord[1]))

    coord = knot2coord(sorted_knots, io_path[-1][1])
    f.write('{0:.2f} {1:.2f}\n'.format(coord[0], coord[1]))

    f.close()


def list_entities(dxf):
    dxf_summary = [shape.dxftype for shape in dxf.entities]
    print('{0:<10}: {1}'.format('LINES', dxf_summary.count('LINE')))
    print('{0:<10}: {1}'.format('ARCS', dxf_summary.count('ARC')))


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


def knots_rank_list_summary(knots_rank):
    print('{0:<16}: {1}'.format('IO knots', [x[1] for x in knots_rank].count(1),
                                [x[0] for x in knots_rank if x[1] == 1]))
    print('{0:<16}: {1}'.format('master knots', [x[1] for x in knots_rank].count(3),
                                [x[0] for x in knots_rank if x[1] == 3]))
    print('{0:<16}: {1}'.format('chain knots', [x[1] for x in knots_rank].count(2),
                                [x[0] for x in knots_rank if x[1] == 2]))


def paths_summary(io_path, ct_path):
    print("-----IN PATHS-----")
    for var in io_path:
        print("{0} {1}".format(var[0], var[1]))

    print("-----CT PATHS-----")
    for var in ct_path:
        print("{0} {1}".format(var[0], var[1]))

    print("-----OUT PATHS-----")
    for var in reversed(io_path):
        print("{0} {1}".format(var[1], var[0]))


def find_path(crit, el_kt_list, sorted_knots, excl_knot):
    path = []
    bar_length = 50
    buff_length = len(el_kt_list)
    sect_length = int(2 * buff_length / bar_length)
    counter = 0
    knots_rank = knots_rank_list(el_kt_list, sorted_knots, excl_knot)
    val_max = max(knots_rank, key=lambda tup: tup[1])[1]
    while val_max > crit:
        knots_rank = knots_rank_list(el_kt_list, sorted_knots, excl_knot)
        curr_knot = knots_rank_find(knots_rank, 1)
        val_max = max(knots_rank, key=lambda tup: tup[1])[1]
        for element in el_kt_list:
            if curr_knot[0] in element:
                counter += 1
         #       print 'found element', element, 'val max', val_max
                if element[0] == curr_knot[0]:
                    path.append(element)
                if element[1] == curr_knot[0]:
                    path.append(element[::-1])

                el_kt_list.remove(element)
                break
            if counter > sect_length:
                counter = 0
    if crit == 1:
        path.append([path[-1][1], path[0][0]])
    return path


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
    if cross_prod(u, v) > 0:
        return False
    else:
        return True


def print_list(data_list, common_text):
    for var in data_list:
        print('{0} {1}'.format(common_text, var))

def interp(s1, s2, p):
    s1x, s1y = s1[:-1]
    s2x, s2y = s2[:-1]
    return [p, (s2y - s1y) / (s2x - s1x) * (p - s1x) + s1y, 0]

def interp_series(x, y, p):
    if x[0] >= p:
        return interp([x[0], y[0], 0], [x[1], y[1], 0], p)[1]
    if x[-1] < p:
        return interp([x[-2], y[-2], 0], [x[-1], y[-1], 0], p)[1]
    for i in range(len(x) - 1):
        if x[i] < p <= x[i + 1]:
            return interp([x[i], y[i], 0], [x[i + 1], y[i + 1], 0], p)[1]

def cross(u, v, opt):
    i = u[1] * v[2] - u[2] * v[1]
    j = u[2] * v[0] - u[0] * v[2]
    k = u[0] * v[1] - u[1] * v[0]
    if opt in 'normalize':
        vec_len = (i**2 + j**2 + k**2)**0.5
        if vec_len > 0:
            i /= vec_len
            j /= vec_len
            k /= vec_len
    return [i, j, k]

def vect(p1, p2):
    u1 = p2[0] - p1[0]
    u2 = p2[1] - p1[1]
    u3 = p2[2] - p1[2]
    return [u1, u2, u3]


def knots2face(knots, t, offset):
    knots_list = knots[:]

    if knots_list[0] != knots_list[-1]:
        knots_list.append(knots_list[0])

    points_list = [Base.Vector(var[0], var[1], var[2]) for var in knots_list]
    poly = Part.makePolygon(points_list)
    face = Part.Face(poly)

    if t != 0:
        n_vect=[0,0,0]
        shift=0
        max_shift=len(knots_list)
        # for var in knots_list:
        #     print var

        while sum(n_vect) ==0:
            n_vect = cross(vect(knots_list[shift], knots_list[1+shift]),
                           vect(knots_list[shift], knots_list[-1+shift]), 'normalize')
            shift+=1
        #    print shift, n_vect, sum(n_vect), vect(knots_list[0], knots_list[1+shift]), vect(knots_list[0], knots_list[-1-shift])

        face = face.extrude(Base.Vector(
            n_vect[0] * -t, n_vect[1] * -t, n_vect[2] * -t))

    return face


def interp_series(x, y, p):
    if x[0] >= p:
        return interp([x[0], y[0], 0], [x[1], y[1], 0], p)[1]
    if x[-1] < p:
        return interp([x[-2], y[-2], 0], [x[-1], y[-1], 0], p)[1]
    for i in range(len(x) - 1):
        if x[i] < p <= x[i + 1]:
            return interp([x[i], y[i], 0], [x[i + 1], y[i + 1], 0], p)[1]


def dist(x1, y1, x2, y2):
    return ((x2 - x1)**2 + (y2 - y1)**2)**0.5


def interp_len(x, y):  # converts x coord to a path distance from starting point s
    path_points = zip(x, y)
    path_segments_point = zip(path_points, path_points[1::])
    path_segments_lengt = [dist(s1[0], s1[1], s2[0], s2[1])
                           for s1, s2 in path_segments_point]

    path_dist = 0.0
    path_segments = [0]
    for seg_len in path_segments_lengt:
        path_dist += seg_len
        path_segments.append(path_dist)

    return path_segments


def balance_knots(prof_knots, fold_knots):
    pool_list = []
    #=============OK section
    for s1, s2 in zip(prof_knots, prof_knots[1::]):
        fit_candidates = [var[0] for var in fold_knots if s1[
            0] < var[0] < s2[0] or s1[0] > var[0] > s2[0]]
        print s1, s2
        pool_list.append(s1)
        if fit_candidates:
            if s2[0] >= s1[0]:
                for extra_point in fit_candidates:
                    point = interp(s1, s2, extra_point)
                    pool_list.append(point)
            else:
                for extra_point in fit_candidates[::-1]:
                    point = interp(s2, s1, extra_point)
                    pool_list.append(point)

    pool_list.append(prof_knots[-1])
    print '\ngenerating pool_list knots: DONE'
    return pool_list
#*********************************************************************DEFAULT PARAMETERS
dflt_o_f = 'all'  # decimal accuracy
dflt_dec_acc = 3  # decimal accuracy
dflt_n_arc = 10  # number of segments
dflt_l_arc = 1  # minimal segment length
dflt_path_dir = 1  # closed path collecting direction
dflt_fc_model_thickness = 0  # default sheet thickness
#*********************************************************************PROGRAM

parser = argparse.ArgumentParser(description='knots lists manipulation kit')
parser.add_argument('-i1', '--master_input_knots', type=str,
                    required=True, help='base profile')
parser.add_argument('-i2', '--slave_input_knots', type=str,
                    help='folding profile')

parser.add_argument('-o', '--output_knots', type=str, help='output file name')
parser.add_argument('-fc', '--fc_model_output_name',
                    type=str, help='generate FreeCAD model')
parser.add_argument('-fct', '--fc_model_thickness', default=dflt_fc_model_thickness,
                    type=float, help='FreeCAD model sheet thickness')

# FreeCAD creates a group UNFOLD with and an unfolded geometry
args = parser.parse_args()
fc_model_name = args.fc_model_output_name
sheet_thickness = args.fc_model_thickness
prof_f = args.master_input_knots #base profile
fold_f = args.slave_input_knots  #folding profile
output_f = args.output_knots

if not output_f:
    output_f = ''.join(['mod_', prof_f])

prof_knots = read_data(prof_f, True)

if not prof_knots:
    print('no knots found, exit program')

else:

    prof_Xmin = min([var[0] for var in prof_knots])
    prof_Xmax = max([var[0] for var in prof_knots])
    prof_Ymin = min([var[1] for var in prof_knots])
    prof_Ymax = max([var[1] for var in prof_knots])
    print('base profile - {0} |X: {1} <-> {2} |Y: {3} <-> {4} '.format(prof_f,prof_Xmin, prof_Xmax,prof_Ymin, prof_Ymax))

    if fold_f:
        fold_knots = read_data(fold_f, True)
    else:
        fold_knots = [[prof_Xmin, 0, 0], [prof_Xmax, 0, 0]]

    X = [var[0] for var in fold_knots]
    Y = [var[1] for var in fold_knots]

    pool_list = balance_knots(prof_knots, fold_knots)

    x_len = interp_len(X, Y)
    x0_len = interp_series(X, x_len, pool_list[0][0])

    unfold_faces_list = [
        [interp_series(X, x_len, p1[0]) - x0_len, p1[1], 0] for p1 in pool_list]

    write_data(output_f, unfold_faces_list, True)

    if fc_model_name:

        print 'generating FreeCAD model: ', fc_model_name

        fc_faces_list = []
        fc_fold_faces_list = []
        fc_unfold_faces_list = []
        folds_X = []

        folds_list = [
            var for var in fold_knots if prof_Xmin < var[0] < prof_Xmax]

        folds_X.append(prof_Xmin)
        folds_X.extend([var[0] for var in folds_list])
        folds_X.append(prof_Xmax)

        for b1, b2 in zip(folds_X, folds_X[1::]):
            fc_faces_list.append(
                [var for var in pool_list if b1 <= var[0] <= b2])

        for face in fc_faces_list:
            fc_fold_faces_list.append(
                [[p1[0], p1[1], interp_series(X, Y, p1[0])] for p1 in face])
            fc_unfold_faces_list.append(
                [[interp_series(X, x_len, p1[0]) - x0_len, p1[1], 0] for p1 in face])

        doc = FreeCAD.newDocument(fc_model_name)
        myPart = doc.addObject('Part::Feature', 'fold_sheet')
        myPart.Shape = Part.Compound(
            [knots2face(face, sheet_thickness, 0) for face in fc_fold_faces_list])
        # doc.saveAs(fc_model_name+'_folded.fcstd')

        myPart = doc.addObject('Part::Feature', 'u_fold_sheet')
        myPart.Shape = Part.Compound(
            [knots2face(face, sheet_thickness, 0) for face in fc_unfold_faces_list])
        doc.saveAs(fc_model_name + '.fcstd')
