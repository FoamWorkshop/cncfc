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
        9 - save ct_path, io_path and reversed(io_path) to the output files

        improvements:
            #. fix dxf paths

            #. layer naming convention
                number with a tag in '()':
                    xxxx()
                    xxxx.()
                    xxxx#y()
                    xxxx#y.()

            #. merging dxf drawings
                *layers with the same number before # are merged so:
                    xxxx#y0('path with props 0')
                    xxxx#y1('path with props 1')
                    xxxx#y2('path with props 2')

            #. drawing keywords:
                heating = 10(W)
                angle = 90(deg)
                radius = 100(mm)
                cut_speed = 200(mm/min)
                start/circle center - indicates begining of an unlooped path

            #. data structure:
                r - radius
                a - angle
                z - axial
                s - slope
                c - cutting speed
                p - heating

        drawing options:
        continous lines
        loops
        '''

import os
import argparse
import sys
import dxfgrabber
import numpy as np
import pickle
import cncfclib


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

def ct_len(io_path, sorted_knots):
    coord_list = []
    for var in io_path:
        coord = knot2coord(sorted_knots, var[0])
        coord_list.append((coord[0], coord[1]))

    coord = knot2coord(sorted_knots, io_path[-1][1])
    coord_list.append((coord[0], coord[1]))

    coord_arr=np.array(coord_list)
    l_arr = np.linalg.norm(np.diff(coord_arr,axis=0), axis=1)
    return np.sum(l_arr)

def ct_len_1(section_list):
    coord_arr=np.array(section_list)
    l_arr = np.linalg.norm(np.diff(coord_arr,axis=0), axis=1)
    return np.sum(l_arr)


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
    if np.linalg.norm(np.cross(u, v)) > 0:
        return False
    else:
        return True

def main(args):

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

            if len(layer_list):
                layer_name_list= [ var.name for var in dxf_layers if var.name in layer_list]
            else:
                layer_name_list = [var.name for var in dxf_layers if not ('~' in var.name or len(var.name)==1)]

            layer_name_list= [ var.name for var in dxf_layers if layer_list[0] in var.name]

            for layer_name in sorted(layer_name_list):
                knots_list, elements_list, segment_bounds, shape_count, start_coord = cncfclib.dxf_read(dxf, layer_name, dec_acc, n_arc, l_arc)
                print 'dxf loaded'

                sorted_knots = knots_dict(knots_list)
                el_kt_list = elements_coords2knots(elements_list, sorted_knots)
                knots_rank = knots_rank_list(el_kt_list, sorted_knots, None)
                master_knot = knots_rank_find(knots_rank, 3)
                IO_knot = knots_rank_find(knots_rank, 1)
                print(master_knot)
                if len(start_coord) and len(IO_knot) % 2 == 0 and master_knot[0] is None:
                    print('found {} lines'.format(len(IO_knot)//2))


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
                                section = z.tofileslist()
                            updated_section_list.append(section)
    #flatten the list of lists
                        section_list = [var for sublist in updated_section_list for var in sublist]
    #SUMMARY
                    speed = 60/200
                    print('{0:11}: {1:4d} | cut len: {2:4.0f} | cut time {3:4.0f}s'.format('i/o  seg.', len(section_list)-1, ct_len_1(section_list), ct_len_1(section_list)*speed))
                    print('{0:11}: {1:4d} | cut len: {2:4.0f} | cut time {3:4.0f}s'.format('loop seg.', len(ct_path),        ct_len(ct_path, sorted_knots), ct_len(ct_path, sorted_knots)*speed))
                    print('{0}'.format('-' * 80))
    #SUMMARY
                    if '1' in output_path:
                        i_file_name = '{1}{2}.{3}'.format(case_name[0], layer_name, '1', 'knt')
                        knots2file_1(i_file_name, section_list, z_coord)
    #                     size = len(section_list)
    #                     a_arr=np.zeros((1,size,1))
    #                     r_arr=np.zeros((1,size,1))
    #                     z_arr=np.zeros((1,size,1))
    #                     v_arr=np.zeros((1,size,2))
    #                     print('found: ', len(section_list),' input sections')
    #                     for i, var in enumerate(section_list):
    #                         a_arr[0,i,0]=var[0]
    #                         r_arr[0,i,0]=var[1]
    #                         z_arr[0,i,0]=0
    #                         v_arr[0,i,:]=np.array([1,0])
    #
    #                     res_dict = {'a_arr':np.rot90(a_arr, k=-1), #rotation angle
    #                                 'r_arr':np.rot90(r_arr, k=-1), #radius R/X
    #                                 'z_arr':np.rot90(z_arr, k=-1), #height Z/Y
    #                                 'v_arr':np.rot90(v_arr, k=-1)} #slope (useful for tapered wings)
    #
    #                     with open(i_file_name, 'wb') as f:
    # # Pickle the 'data' dictionary using the highest protocol available.
    #                         pickle.dump(res_dict, f, pickle.HIGHEST_PROTOCOL)

                    if '3' in output_path:
                        o_file_name = '{1}{2}.{3}'.format(case_name[0], layer_name, '3', 'knt')
                        knots2file_1(o_file_name, section_list[::-1], z_coord)
                        # size = len(section_list)
                        # a_arr=np.zeros((1,size,1))
                        # r_arr=np.zeros((1,size,1))
                        # z_arr=np.zeros((1,size,1))
                        # v_arr=np.zeros((1,size,2))
                        # print('found: ', len(section_list),' output sections')
                        # for i, var in enumerate(section_list):
                        #     a_arr[0,i,0]=var[0]
                        #     r_arr[0,i,0]=var[1]
                        #     z_arr[0,i,0]=0
                        #     v_arr[0,i,:]=np.array([1,0])
                        #
                        # res_dict = {'a_arr':np.rot90(a_arr, k=-1), #rotation angle
                        #             'r_arr':np.rot90(r_arr, k=-1), #radius R/X
                        #             'z_arr':np.rot90(z_arr, k=-1), #height Z/Y
                        #             'v_arr':np.rot90(v_arr, k=-1)} #slope (useful for tapered wings)
                        #
                        # with open(o_file_name, 'wb') as f:
                        # # Pickle the 'data' dictionary using the highest protocol available.
                        #     pickle.dump(res_dict, f, pickle.HIGHEST_PROTOCOL)

                    if '2' in output_path:
                        ct_file_name = '{1}{2}.{3}'.format(case_name[0], layer_name, '2', 'knt')
                        knots2file(ct_file_name, ct_path, sorted_knots)
                        # size = len(section_list)+1
                        # a_arr=np.zeros((1,size,1))
                        # r_arr=np.zeros((1,size,1))
                        # z_arr=np.zeros((1,size,1))
                        # v_arr=np.zeros((1,size,2))
                        #
                        # print('found: ', len(section_list)+1,' shape sections')
                        # for i, var in enumerate(ct_path):
                        #     coord = knot2coord(sorted_knots, var[0])
                        #     # print(coord)
                        #     a_arr[0,i,0] = 0
                        #     r_arr[0,i,0] = coord[0]
                        #     z_arr[0,i,0] = coord[1]
                        #     v_arr[0,i,:] = np.array([1,0])
                        #
                        # a_arr[0,-1,0]=a_arr[0,0,0]
                        # r_arr[0,-1,0]=r_arr[0,0,0]
                        # z_arr[0,-1,0]=0
                        # v_arr[0,-1,:]=v_arr[0,0,:]
                        #
                        # res_dict = {'a_arr':np.rot90(a_arr, k=-1), #rotation angle
                        #             'r_arr':np.rot90(r_arr, k=-1), #radius R/X
                        #             'z_arr':np.rot90(z_arr, k=-1), #height Z/Y
                        #             'v_arr':np.rot90(v_arr, k=-1)} #slope (useful for tapered wings)
                        #
                        # with open(ct_file_name, 'wb') as f:
                        # # Pickle the 'data' dictionary using the highest protocol available.
                        #     pickle.dump(res_dict, f, pickle.HIGHEST_PROTOCOL)

                    print(' saved')

if __name__ == '__main__':
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
    parser.add_argument('-larc', '--arc_seg_len', type=float,
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

    parser.add_argument('-plt', '--plot_paths', action='store_true', help='plot cutting paths')

    args = parser.parse_args()



    # parser = argparse.ArgumentParser(description='This program slices dxf model using XYZ planes and saves projection of the section curves to a dxf file. Each section curve is saved as a block inserted into corresponding layer. In order to make the usage easier, block are colored.')
    # parser.add_argument('-i', '--input', required = True, type=str, help='stl input file')
    # parser.add_argument('-o', '--output', default=[], type=str, help='dxf output file')
    # parser.add_argument('-s', '--slices', default=[10,1,1], nargs='+', type=int, help='(3 int numbers) number of slices per axis.\n 1 section -> mid plane\n2 - end section planes, >2 - linear spacing between min and max dimension')
    # parser.add_argument('-k_e', '--slices_end_factors', default=[1,1,1], nargs='+', type=float, help='end sections offset (3 float numbers)')
    # parser.add_argument('-p_o', '--proj_offsets', action='store_true', help='projection offsets. if activated additional coord ')
    # parser.add_argument('-p_s', '--plot_sections', action='store_true', help='plot 3d sections using matplotlib')
    # parser.add_argument('-sx', '--x_slices', default=[], nargs='+', type=float, help='define x slices, if -s is specified, the 2 lists are merged')
    # parser.add_argument('-sy', '--y_slices', default=[], nargs='+', type=float, help='define y slices, if -s is specified, the 2 lists are merged')
    # parser.add_argument('-sz', '--z_slices', default=[], nargs='+', type=float, help='define z slices, if -s is specified, the 2 lists are merged')
    # args = parser.parse_args()
    # print(args)
    # fname_stl = args.input
    # fname_dxf = args.output
    # slices_list = args.slices
    # k_e = args.slices_end_factors
    # proj_offsets = args.proj_offsets
    # plot_sections = args.plot_sections
    main(args)
