#!/usr/bin/python3
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
            #. use regular expresion to match layer names
            #. layer naming convention
                number with a tag in '()':
                    xxxx()
                    xxxx.()
                    xxxx#y()
                    xxxx#y.()

            #. merging dxf drawings
                *layers with the same number before # are merged so:
                    xxxx#a('path with props 0')
                    xxxx#b('path with props 1')
                    xxxx#1('path with props 2')

            #. makeing drafting profiles
                *layers with the same number before # are merged so:
                    xxxx$0('path with props 0')
                    xxxx$1('path with props 1')

            #. drawing keywords:
                heating = 10(W)
                angle = 90(deg)
                radius = 100(mm)
                feed = 200(mm/min)
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
        import os
        '''

import argparse
import sys
import dxfgrabber
import numpy as np
import pickle
import cncfclib
import os
import re
import matplotlib.pyplot as plt

def ct_sect(sect_arr):
    u = sect_arr[:,0,:]
    v = sect_arr[:,1,:]
    p = v - u
    l_arr = np.linalg.norm( p, axis=1)
    return l_arr

def ct_len_1(l_arr):
    return np.sum(ct_sect(l_arr))

def ct_speed(l_arr, prop, prop_dict):
    sect_speed = ct_sect(l_arr) * 60 / np.array([prop_dict[key]['feed'] for key in prop])
    return np.sum(sect_speed)

def print_stats(ss):
    ct_len_list = []
    ct_time_list = []
    print('{:-^79}'.format('CUTTING STATS'))
    for var1 in ss:
        if var1:
            io_path1, lo_path1, io_prop1, lo_prop1, prop_dict1 = var1[0]
            io_path2, lo_path2, io_prop2, lo_prop2, prop_dict2 = var1[1]
            # print(prop_dict1)
            if prop_dict1[0]['radius'][2]> prop_dict2[0]['radius'][2]:
                master_io_path = io_path1
                master_lo_path = lo_path1
                master_io_prop = io_prop1
                master_lo_prop = lo_prop1
                master_prop_dict = prop_dict1
            else:
                master_io_path = io_path2
                master_lo_path = lo_path2
                master_io_prop = io_prop1
                master_lo_prop = lo_prop1
                master_prop_dict = prop_dict2

            io_ct_len = ct_len_1(master_io_path)
            lo_ct_len = ct_len_1(master_lo_path)
            io_ct_speed = ct_speed(master_io_path, master_io_prop, master_prop_dict)
            lo_ct_speed = ct_speed(master_lo_path, master_lo_prop, master_prop_dict)

            print('layer:  {}'.format(master_prop_dict[0]['layer']))
            # print('radius:  {}'.format(master_prop_dict[0]['radius']))
            print('io cutting len:  {0:10.0f}mm'.format(io_ct_len))
            print('lo cutting len:  {0:10.0f}mm'.format(lo_ct_len))
            print('io cutting time: {0:10.0f}s'.format(io_ct_speed))
            print('lo cutting time: {0:10.0f}s'.format(lo_ct_speed))
            ct_len_list.append([io_ct_len, lo_ct_len])
            ct_time_list.append([io_ct_speed+lo_ct_speed])

    total_ct_len = np.sum(np.array(ct_len_list))
    total_ct_time = np.sum(np.array(ct_time_list))
    print('{:-^79}'.format('CUTTING SUMMARY'))
    print('total cutting length: {:5.0f}mm'.format(total_ct_len))
    print('total cutting time:   {:5.0f}s'.format(total_ct_time))
def make_path_chain(io,lo):
    return np.vstack((io[:,0,:],    io[-1,1,:],
                      lo[:,0,:],    lo[-1,1,:],
                      io[::-1,1,:], io[0,0,:]))

def make_prop_chain(io,lo):
    return np.hstack((io,    io[-1],
                      lo,    lo[-1],
                      io[::-1], io[0]))
def make_projection(chain1, chain2):
    return chain1[:,:2], chain2[:,:2]

def cut_projection2gcode(cut_path1, cut_path2, cut_prop1, cut_prop2, prop_dict1, prop_dict2, machine_conf, gcode_conf):
    gcode=[]
    speed_old = 0
    for cp1, cp2, c_prop1 in zip(cut_path1, cut_path2, cut_prop1):
        speed_new=prop_dict1[c_prop1]['feed']
        angle = prop_dict1[c_prop1]['angle']
        if speed_new != speed_old:
            line = '{0[0]} {1}{0[1]}'.format(gcode_conf['spindle'],
                                    speed_new)
            gcode.append(line)
            speed_old = speed_new


        line = 'G1 {1[0]} {0[0][0]:<8.2f} ' \
                  '{1[1]} {0[0][1]:<8.2f} ' \
                  '{1[2]} {0[1][0]:<8.2f} ' \
                  '{1[3]} {0[1][1]:<8.2f} ' \
                  '{1[4]} {0[2]:<8.2f} '.format([cp1, cp2, angle], machine_conf['ax'])
        gcode.append(line)

    gcode.append('{0[2]}'.format(gcode_conf['spindle']))

    return gcode
def print_gcode(ss):
    machine_conf={
#axA master column
    'ax':['X','Y', 'U', 'V', 'B'],
#distance between columns.
    'AB_dist':0.45,
# distance between rotary table and column A
    'AT_dist':0.225}
    gcode_conf={'comment':['(',')'],
                'spindle':['M3','S','M5']}

    ct_len_list = []
    ct_time_list = []

    print('{:-^79}'.format('GCODE'))
    for i, var1 in enumerate(ss):
        if var1:
            io_path1, lo_path1, io_prop1, lo_prop1, prop_dict1 = var1[0]
            io_path2, lo_path2, io_prop2, lo_prop2, prop_dict2 = var1[1]

            print('(file :{})'.format(None))
            print('(layers : {0[0]} {0[1]})'.format([prop_dict1['layer'], prop_dict2['layer']]))
            print('{0[0]}sequence :{1}{0[1]}'.format(gcode_conf['comment'],i))


            chain1_path = make_path_chain(io_path1, io_path1)
            chain2_path = make_path_chain(io_path2, io_path2)
            # print(io_prop1)
            # print(lo_prop1)
            chain1_prop = make_prop_chain(io_prop1, lo_prop1)
            chain2_prop = make_prop_chain(io_prop2, lo_prop2)

            cut_path1, cut_path2 = make_projection(chain1_path, chain2_path)

            gcode = cut_projection2gcode(cut_path1, cut_path2, chain1_prop, chain2_prop, prop_dict1, prop_dict2, machine_conf, gcode_conf)
            for line in gcode:
                print(line)


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

    req_layer = layer_list[0]
    for i, files_dxf_member in enumerate(files_dxf):

        case_name = os.path.splitext(files_dxf_member)
        dxf = dxfgrabber.readfile(files_dxf_member, {"assure_3d_coords": True})
        dxf_layers = [var.name for var in dxf.layers]
        regex2 = re.compile("^({})#.*".format(req_layer), re.IGNORECASE)
        z2 = [layer for layer in dxf_layers for m in [regex2.search(layer)] if m]
        dxf_params = (dec_acc, n_arc, l_arc)

        if z2:
            print(z2)
            print('sequences')
            regex3 = re.compile("^\w+#\d+", re.IGNORECASE)
            z3 = sorted(list(set([m.group() for layer in z2 for m in [regex3.search(layer)] if m])))
            print(z3)
            seq_list =[]
            for seq_step in z3:
                regex4 = re.compile("{}#[0-1]".format(seq_step), re.IGNORECASE)
                z5 = sorted(list(set([m.group() for layer in z2 for m in [regex4.search(layer)] if m])))
                plane_sections_list = []
                for z6 in z5:
                    regex5 = re.compile("{}#.*".format(z6), re.IGNORECASE)
                    z7 = sorted(list(set([m.group() for layer in z2 for m in [regex5.search(layer)] if m])))
                    plane_sections_list.append(z7)
                seq_list.append(plane_sections_list)
            ss=[]
            for i, seq in enumerate(seq_list):
                print('\nseq num: ', i)
                pp0 = []
                pp1 = []
                pp =[]
                for j, plane in enumerate(seq):
                    print(' plane: ', j)
                    for section in plane:
                        print('  ',section)

                    io_path1, lo_path1, io_path_prop1, lo_path_prop1, prop_dict1 = cncfclib.extract_dxf_path(dxf, plane, dxf_params)
                    pp.append([io_path1, lo_path1, io_path_prop1, lo_path_prop1, prop_dict1])

                if len(pp)==1:
                    ss.append(pp*2)
                else:
                    ss.append(pp)

        print_stats(ss)
        print_gcode(ss)

            # cncfclib.plot_path1(ss)
            # print(['o']*10)


#SUMMARY
#                 if '1' in output_path:
#                     i_file_name = '{1}{2}.{3}'.format(case_name[0], layer_name, '1', 'knt')
#                     np.save(i_file_name, io_path)
#
#                 if '3' in output_path:
#                     o_file_name = '{1}{2}.{3}'.format(case_name[0], layer_name, '3', 'knt')
#                     np.save(o_file_name, io_path[::-1])
#
#                 if '2' in output_path:
#                     ct_file_name = '{1}{2}.{3}'.format(case_name[0], layer_name, '2', 'knt')
#                     # knots2file(ct_file_name, ct_path, sorted_knots)
#                     # size = len(section_list)+1
#                     # a_arr=np.zeros((1,size,1))
#                     # r_arr=np.zeros((1,size,1))
#                     # z_arr=np.zeros((1,size,1))
#                     # v_arr=np.zeros((1,size,2))
#                     #
#                     # print('found: ', len(section_list)+1,' shape sections')
#                     # for i, var in enumerate(ct_path):
#                     #     coord = knot2coord(sorted_knots, var[0])
#                     #     # print(coord)
#                     #     a_arr[0,i,0] = 0
#                     #     r_arr[0,i,0] = coord[0]
#                     #     z_arr[0,i,0] = coord[1]
#                     #     v_arr[0,i,:] = np.array([1,0])
#                     #
#                     # a_arr[0,-1,0]=a_arr[0,0,0]
#                     # r_arr[0,-1,0]=r_arr[0,0,0]
#                     # z_arr[0,-1,0]=0
#                     # v_arr[0,-1,:]=v_arr[0,0,:]
#                     #
#                     # res_dict = {'a_arr':np.rot90(a_arr, k=-1), #rotation angle
#                     #             'r_arr':np.rot90(r_arr, k=-1), #radius R/X
#                     #             'z_arr':np.rot90(z_arr, k=-1), #height Z/Y
#                     #             'v_arr':np.rot90(v_arr, k=-1)} #slope (useful for tapered wings)
#                     #
#                     # with open(ct_file_name, 'wb') as f:
#                     # # Pickle the 'data' dictionary using the highest protocol available.
#                     #     pickle.dump(res_dict, f, pickle.HIGHEST_PROTOCOL)

    print('\nDone. Thank you!')

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

    main(args)
