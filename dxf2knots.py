#!/usr/bin/python
from __future__ import division
__author__ = 'FoamWorkshop'

'''
program options:
-i [input files | all]
-a [float decimal accuracy] default 3
-narc [number o segments] - default 10
-larc [minimal segment length] - default 1
-cw   [1|0] default 1: 1 - clockwise; 0 - counter clockwise of the  closed path
-l process selected layers

info:
    the program extracts entities:
        LINE, ARC
    and converts them into a list of coordinates defining:
        1 - input/output path
        2 - closed contour paths

path finding alghoritm:
        1 - find a segment including the start knot (k0)
        2 - read the remining knot (r_k)
        3 - push (remove from the pool) the segment knots to the path list ordered like: [k0, r_k]
        4 - assign the new start knot as (k0) and repeat from 1 while the pool is not empty.
        5 - if the closed contour path, add the last segment connecting the first knot and the last in the path list [k(-1), k(0)]

program algorithm:
        1 - open dxf and go through all layers
        2 - find the entities and read DATA with specified accuracy:
             LINE - read START, END coordinates
             ARC - read START ANGLE, END ANGLE and RADIUS:
             convert the curve to segments with following options:
              n - number of segments
              or l_min - specified minimal segment length

            note: the concept is taken from Abaqus inp. files organisation
            the  output is:
             - list of segment coordinates [[x1,y1,y1], [x2,y2,y2];...]
             - list of segment knots [[x1,y1,y1]; [x2,y2,y2]; [xn, yn, zn];...]

        3 - remove duplicates from the segment knots list. position of a coordinate int the knot list indicates the knot number.
        4 - replace the segment list coordinates by knot numbers. the resultant list includes segment knot numbers [[k1, k2];[ki kn];...]
        5 - sort the segment list by knots count. the proper list should include 1 - io knot shared with a segment and: 1 - master knot shared between 3 segments
        6 - find the io_path with begin in io knot and end in master knot. the outcome is:
                io_path and remining segment pool for the closed path
        7 - find the END segment for clock wise or counter clock wise direction and exclude the last knot from the ranking list.
        8 - find the ct_path
        9 - save ct_path, io_path and reversed(io_path) to the output files'''

import os
import argparse
import sys
import dxfgrabber
import numpy as np
from cncfclib import *


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
            knots_rank.append([var[0], [x for a in el_kt_list for x in a].count(var[0])])
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
        f.write('{0:.3f} {1:.3f}\n'.format(coord[0], coord[1]))

    coord = knot2coord(sorted_knots, io_path[-1][1])
    f.write('{0:.3f} {1:.3f}\n'.format(coord[0], coord[1]))

    f.close()

def knots2file_1(name, section_list, z_coord):
    f = open(name, 'w')
    # print('TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT')
    for var in section_list:
        # coord = knot2coord(sorted_knots, var[0])
        # print(var)
        f.write('{0:.3f} {1:.3f} {2:.3f}\n'.format(var[0], var[1], z_coord))

    # coord = knot2coord(sorted_knots, io_path[-1][1])
    # f.write('{0:.3f} {1:.3f}\n'.format(coord[0], coord[1]))

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

    knots_rank = knots_rank_list(el_kt_list, sorted_knots, excl_knot)

    curr_knot = knots_rank_find(knots_rank, 1)
    last_knot = knots_rank_find(knots_rank, 3)

    curr_element=[]

    # print '\nfinding path'
    while not ((curr_element is None) or curr_knot[0]==last_knot[0]):
        # print '\rpool size: {0}'.format(len(path)),

        curr_element=next((element for element in el_kt_list if curr_knot[0] in element), None)
        if not (curr_element is None):

            if curr_element[0] == curr_knot[0]:
                curr_knot=[curr_element[1]]
                path.append(curr_element)
            else:
                curr_knot=[curr_element[0]]
                path.append(curr_element[::-1])

            el_kt_list.remove(curr_element)

    if crit == 1:
        path.append([path[-1][1], path[0][0]])
    # print '\n'
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


def dxf_read(files, layer_name, dec_acc, n_arc, l_arc):
    tol = dec_acc
    knots_list = []
    elements_list = []
    hrd_knots_list=[]
    hrd_element_list=[]
    segment_bounds=[]
    line_count = 0
    arc_count = 0
    circle_count = 0
    path_offset =[0,0,0]
#    list_entities(dxf)

    for shape in dxf.entities:
        # print(shape.dxftype)
        if shape.layer == layer_name:
            if shape.dxftype == 'SPLINE':

                print('degree: ',shape.degree)
                print('start tangent: ',shape.start_tangent)
                print('end tangent: ',shape.end_tangent)
                print('control points: ',shape.control_points)
                print('fit points: ',shape.fit_points)
                print('knots: ',shape.knots)
                print('weights: ',shape.weights)
                print('normal vector: ',shape.normal_vector)

                # circle_count += 1
                # p1 = tuple(round(x, tol) for x in shape.center)
                # segment_bounds.append(p1)
                # # print(circle_count, p1)

            if shape.dxftype == 'CIRCLE':
                circle_count += 1
                p1 = tuple(round(x, tol) for x in shape.center)
                segment_bounds.append(p1)
                # print(circle_count, p1)

            if shape.dxftype == 'LINE':
                line_count += 1
                p1 = tuple(round(x, tol) for x in shape.start)
                p2 = tuple(round(x, tol) for x in shape.end)
                if p1!=p2:
                    knots_list.append(p1)
                    knots_list.append(p2)
                    elements_list.append([p1, p2])

            if shape.dxftype == 'MTEXT':
                if shape.raw_text == 'coord_0':
                    # print('path offset: {}'.format(shape.insert))
                    path_offset = tuple(round(x, tol) for x in shape.insert)
                    print('path offset: {}'.format(path_offset))

            if shape.dxftype == 'ARC':
                arc_count += 1
                ARC_knots_list = []
                n = n_arc  # number of segments
                min_len = l_arc
                O = shape.center
                R = shape.radius
                angl_1 = shape.start_angle * np.pi / 180
                angl_2 = shape.end_angle * np.pi / 180

                if angl_2 >= angl_1:
                    angl_list = np.linspace(angl_1, angl_2, n)
                else:
                    angl_list = np.linspace(angl_1, angl_2 + 2 * np.pi, n)

                arc_len = R * np.absolute(angl_2 - angl_1)

                if arc_len / n < min_len:
                    n = max(int(arc_len / min_len), 3)

                for angl in angl_list:
                    ARC_knots_list.append(
                        (round(O[0] + R * np.cos(angl), tol), round(O[1] + R * np.sin(angl), tol), O[2]))

                for i in range(n - 1):
                    elements_list.append(ARC_knots_list[i:i + 2])

                knots_list.extend(ARC_knots_list)

    # print 'remove duplicates'
    # print 'number of segments: ', len(elements_list)

    hrd_element_list.append(elements_list[0])

    for var in elements_list:
        tmp=[var for hrd_var in hrd_element_list if (var[0] in hrd_var) and (var[1] in hrd_var)]
        if  not len(tmp):
            hrd_element_list.append(var)
            # print 'removed elemts: ', len(elements_list)-len(hrd_element_list)


    knots_list=[(var[0]-path_offset[0], var[1]-path_offset[1], var[2]-path_offset[2]) for var in knots_list]
    hrd_element_list=[[(var1[0]-path_offset[0], var1[1]-path_offset[1], var1[2]-path_offset[2]), (var2[0]-path_offset[0], var2[1]-path_offset[1], var2[2]-path_offset[2])]for var1, var2 in hrd_element_list]
    segment_bounds = [(var[0]-path_offset[0], var[1]-path_offset[1], var[2]-path_offset[2]) for var in segment_bounds]
    return (knots_list, hrd_element_list, segment_bounds, [line_count, arc_count])


#*********************************************************************DEFAULT PARAMETERS
dflt_dxf_list = 'all'
dflt_dec_acc = 4  # decimal accuracy
dflt_n_arc = 10  # number of segments
dflt_l_arc = 0.1  # minimal segment length
dflt_path_dir = 1  # closed path collecting direction
#*********************************************************************PROGRAM
parser = argparse.ArgumentParser(description='test')
parser.add_argument('-i', '--input', nargs='+', help='input filenames')
parser.add_argument('-l', '--layer', nargs='+', help='input layers')
parser.add_argument('-a', '--accuracy', type=int,
                    default=dflt_dec_acc, help='decimal accuracy, default: 3')
parser.add_argument('-narc', '--arc_seg_num', type=int,
                    default=dflt_n_arc, help='arc segments number, default: 10')
parser.add_argument('-larc', '--arc_seg_len', type=int,
                    default=dflt_l_arc, help='minimal arc segment length, default: 0.1')
parser.add_argument('-cw', '--collection_dir', type=int,
                    default=dflt_path_dir, help='closed path collection dir')
parser.add_argument('-eq', '--equivalence_knots', type=int,
                    default=False, help='equivalence knots sections to specified number')

parser.add_argument('-eq_skip', '--skip_eq_sections', nargs='+', type=int,
                    default=[], help='equivalence knots sections to specified number')

parser.add_argument('-z', '--z_coord', type=float,
                    default=0, help='add z coordinate to knots')

parser.add_argument('-op', '--output_path', type=str,
                    default='123', help='output path request')

args = parser.parse_args()

dxf_list = args.input
layer_list = args.layer
dec_acc = args.accuracy
n_arc = args.arc_seg_num
l_arc = args.arc_seg_len
path_dir = args.collection_dir
eq_sect = args.equivalence_knots
eq_sect_skip = args.skip_eq_sections
z_coord = args.z_coord
output_path = args.output_path

files_dxf = dxf_list

if 1:

    print('SETTINGS:')
    print('{0}{1:<30}: {2}'.format(' ' * 10, 'decimal accuracy', dec_acc))
    print('{0}{1:<30}: {2}'.format(' ' * 10, 'arc segments count', n_arc))
    print('{0}{1:<30}: {2}'.format(' ' * 10, 'minimal arc segment length', l_arc))
    print('{0}{1:<30}: {2}'.format(' ' * 10, 'equivalence sections', eq_sect))
    print('{0}{1:<30}: {2}'.format(' ' * 10, 'skip equivalence sections', eq_sect_skip))
    print('{0}{1:<30}: {2}'.format(' ' * 10, 'closed path collection dir', path_dir))
    print('{0}{1:<30}: {2}'.format(' ' * 10, 'files', files_dxf))
    print('{0}{1:<30}: {2}'.format(' ' * 10, 'layer name', layer_list[0]))
    print('{0}{1:<30}: {2}'.format(' ' * 10, 'output paths', output_path))

    print('{0}'.format('-' * 80))

    for i, files_dxf_member in enumerate(files_dxf):

        case_name = os.path.splitext(files_dxf_member)
        dxf = dxfgrabber.readfile(files_dxf_member, {"assure_3d_coords": True})
        dxf_layers = dxf.layers


        # print(layer_list)
        # for var in dxf_layers:
        #     print(var.name)
        # if len(layer_list):
        #
        #     layer_name_list= [ var.name for var in dxf_layers if var.name in layer_list]
        #
        #
        # else:
        #     layer_name_list = [var.name for var in dxf_layers if not ('~' in var or len(var.name)==1)]
        #
        # print(layer_name_list)


# <<<<<<< HEAD
        if len(layer_list):
# <<<<<<< HEAD
            # layer_name_list= [ var.name for var in dxf_layers if var in layer_list]
# =======
            layer_name_list= [ var.name for var in dxf_layers if var.name in layer_list]
# >>>>>>> f979f5e... added loft cut
        else:
            layer_name_list = [var.name for var in dxf_layers if not ('~' in var.name or len(var.name)==1)]
# =======

# =======
# >>>>>>> 9d1cd0a... added option eq, eq_skip working eqbalance sections
        layer_name_list= [ var.name for var in dxf_layers if layer_list[0] in var.name]
# >>>>>>> ad63914... added coord_0 text as a ref coordinate system

        for layer_name in sorted(layer_name_list):
            knots_list, elements_list, segment_bounds, shape_count = dxf_read(dxf, layer_name, dec_acc, n_arc, l_arc)
            print 'dxf loaded'

            sorted_knots = knots_dict(knots_list)
            el_kt_list = elements_coords2knots(elements_list, sorted_knots)
            knots_rank = knots_rank_list(el_kt_list, sorted_knots, None)
            master_knot = knots_rank_find(knots_rank, 3)
            IO_knot = knots_rank_find(knots_rank, 1)


            if len(IO_knot) != 1 or len(master_knot) != 1 or IO_knot[0] == None or master_knot[0] == None:
                print('{0:^20}|'.format('SKIPPED'))
                for var in IO_knot:
                    print("IO knot error: {0} coord: {1}".format(var,knot2coord(sorted_knots, var)))
                for var in master_knot:
                    print("master knot error: {0} coord: {1}".format(var,knot2coord(sorted_knots, var)))
                var =10
                print("master knot error: {0} coord: {1}".format(var,knot2coord(sorted_knots, var)))

            else:

                io_path = find_path(2, el_kt_list, sorted_knots, None)  # IO path
                last_el, excl_knot = find_l_el(path_dir, el_kt_list, sorted_knots, master_knot[0])
                ct_path = find_path(1, el_kt_list, sorted_knots, excl_knot[0])  # loop path


                io_knots_coord = [knot2coord(sorted_knots,var[0]) for var in io_path]
                io_knots_coord.append(knot2coord(sorted_knots,io_path[-1][1]))

#EQUIVALENCE SECTION
                section_list = io_knots_coord

                if eq_sect:
                    print('drw. splits: {0:4d}'.format(len(segment_bounds)))
                    section_list = []
                    section = []
                    updated_section_list = []
                    n_seg = np.linspace(0,1,eq_sect)
                    # print(segment_bounds)
                    for i, var in enumerate(io_knots_coord):
                        section.append(var)
                        # print(var)
                        if var in segment_bounds:
                            # print(var)
                            section_list.append(section)
                            section = []
                            section.append(var)
                    section_list.append(section)

                    for i, section in enumerate(section_list):
                        if i not in eq_sect_skip and i not in [var+len(segment_bounds) for var in eq_sect_skip]:
                            # print('equivalence section {}'.format(i))
                            p = np.array(section)
                            l = np.sqrt(np.diff(p[:,0])**2 + np.diff(p[:,1])**2)
                            l = np.cumsum(np.hstack((0,l)))
                            # print(n_seg)
                            # print(l)
                            l_norm = l/l[-1]
                            # print(l_norm)
                            x=np.interp(n_seg, l_norm, p[:,0])
                            y=np.interp(n_seg, l_norm, p[:,1])
                            z=np.vstack((x,y)).T
                            # print(z)
                            section = z.tolist()
                        updated_section_list.append(section)
                    #flatten the list of lists
                    section_list = [var for sublist in updated_section_list for var in sublist]
#EQUIVALENCE SECTION

#SUMMARY
                print('{0:11}: {1:4d}\n'.format('i/o  seg.', len(section_list)-1)),
                print('{0:11}: {1:4d}\n'.format('loop seg.', len(ct_path)))
                print('{0}'.format('-' * 80))
#SUMMARY
                if '1' in output_path:
                    i_file_name = '{1}{2}.{3}'.format(case_name[0], layer_name, '1', 'knt')
                    knots2file_1(i_file_name, section_list, z_coord)
                # knots2file(i_file_name, io_path, sorted_knots)
                if '3' in output_path:
                    o_file_name = '{1}{2}.{3}'.format(case_name[0], layer_name, '3', 'knt')
                    knots2file_1(o_file_name, section_list[::-1], z_coord)
                # knots2file(o_file_name, [var[::-1]
                #                          for var in io_path[::-1]], sorted_knots)
                if '2' in output_path:
                    ct_file_name = '{1}{2}.{3}'.format(case_name[0], layer_name, '2', 'knt')
                    knots2file(ct_file_name, ct_path, sorted_knots)
