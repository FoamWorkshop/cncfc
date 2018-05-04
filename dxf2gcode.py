#!/usr/bin/python3
from __future__ import division
__author__ = 'FoamWorkshop'

import argparse
import sys
import cncfclib
import gcodelib
from cncfc_obj import chain, AxisProfile, ModelProfile, CuttingSpace

# def print_setings(args):
#
#     dxf_list = args.input
#     layer_list = args.layer
#     dec_acc = args.accuracy
#     n_arc = args.arc_seg_num
#     l_arc = args.arc_seg_len
#     path_dir = args.collection_dir
#     eq_sect = args.equivalence_knots
#     eq_sect_skip = args.skip_eq_sections
#     z_coord = args.z_coord
#     output_path = args.output_path
#
#     print('SETTINGS:')
#     print('{0}{1:<30}: {2}'.format(' ' * 10, 'decimal accuracy', dec_acc))
#     print('{0}{1:<30}: {2}'.format(' ' * 10, 'arc segments count', n_arc))
#     print('{0}{1:<30}: {2}'.format(' ' * 10, 'minimal arc segment length', l_arc))
#     print('{0}{1:<30}: {2}'.format(' ' * 10, 'equivalence sections', eq_sect))
#     print('{0}{1:<30}: {2}'.format(' ' * 10, 'skip equivalence sections', eq_sect_skip))
#     print('{0}{1:<30}: {2}'.format(' ' * 10, 'closed path collection dir', path_dir))
#     print('{0}{1:<30}: {2}'.format(' ' * 10, 'files', files_dxf))
#     print('{0}{1:<30}: {2}'.format(' ' * 10, 'layer name', layer_list[0]))
#     print('{0}{1:<30}: {2}'.format(' ' * 10, 'output paths', output_path))
#     print('{0}'.format('-' * 80))
#     return 0

def main(args):
    """
    DXF model organisation:
        layers
        polyline
        profile
        sections
        cutting sections

    machine oconfiguration:
        *------------>XY
        |    |
        |    | d
        |    |
        |s   O--- 90deg
        |    |
        |    0deg
        |
        *------------>UV

        O - rotary table axis, CCW dir, 0 towards UV axis
        d - rotary table distance from XY
        s - distance  between XY and UV columns

    default configuration:
        {'s': 480, 'd': 240, 'max_feed'}

    layer naming convention:
        A#XX#Y#ZZ~comment
        A - layername
        XX- sequence_number
        Y - column number 0 - XY, 1-UV
        ZZ- section number
        ~comment

    layer options specified int the textbox:
        local (layer level):
            feed - feed rate (NOT SCALED)
            power - power on the wire
            angle - rotary table rotation
            radius - layer distance from the RT axis
            coord_0 - reference coordinate sys, represents RT axis
       TODO cut_dir - cutting dir for closed polygons
       TODO split - split ploly line to number of sections

        global (profile level):
            start - indicates the starting point of the cut
    """

    fname_dxf = args.input
    lname_dxf = args.layer

    seq_list, seq_dict = cncfclib.layers2seq(fname_dxf, lname_dxf)

    conf = {'Z_span': 480, 'RTable_loc': 240}

    ct = CuttingSpace(conf)
    md = ModelProfile()

    for j,  seq_name in enumerate(seq_list):
        ap = AxisProfile()

        for i, axis_name in enumerate(seq_name):
            s1 = chain(fname_dxf)

            for prof_name in axis_name:
                s1.AddSeg(prof_name)
                s1.MakeSplit()

            s1.ApplyTransformations()
            s1.MakeChain()
            # s1.PlotChain()

            ap.Add2Axis(i, s1)
            # ap.Plot(mode='2D')
        md.Add2Prof(j, ap)

    # md.Plot(mode = '3D')
    ct.Add2Cut(md)
    # ct.Plot(mode='3D')
    ct.SaveGcode('test.ngc')

if __name__ == '__main__':
    #*********************************************************************DEFAULT PARAMETERS
    dflt_dxf_list = 'all'
    dflt_dec_acc = 4  # decimal accuracy
    dflt_n_arc = 10   # number of segments
    dflt_l_arc = 0.1  # minimal segment length
    dflt_path_dir = 1 # closed path collecting direction
    #*********************************************************************PROGRAM

    parser = argparse.ArgumentParser(description='test')
    parser.add_argument('-i', '--input', type=str, help='input filenames')
    parser.add_argument('-l', '--layer', type=str, help='input layers')
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
