#!/usr/bin/python3
from __future__ import division
__author__ = 'FoamWorkshop'

import argparse
import configparser
import sys
import cncfclib
import gcodelib
from cncfc_obj import chain, AxisProfile, ModelProfile, CuttingSpace, layers2seq
import re

def main(args):
    """
    (current developement)
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

        XY - axis 0
        UV - axis 1
        O - rotary table axis, CCW dir, 0 towards UV axis
        d - rotary table distance from XY
        s - distance  between XY and UV columns

    the default machine configuration is as follows:
        s=480
        d=240
        max_feed=200

    Since the cutting model is organized by layer names it is important to follow layer naming convention:
        A#XX#Y#ZZ~comment
        A - layername
        XX- sequence_number
        Y - column number 0 - XY, 1-UV
        ZZ- section number
        ~comment - optional parameter with description

    each layer can have separate LAYER options can be changed by adding MTEXT field to the layer with an appropriate assignment:
        option_name_1 = value
        option_name_2 = value
        option_name_n = value

    where avaliable LAYER options are:
            feed[float <0, 100>] - feed rate (NOT SCALED)
            power[float <0, 100>] - power on the wire
            angle[float] - rotary table rotation
            radius[float] - layer distance from the RT axis
            coord_0[float, float, float] - reference coordinate sys, represents RT axis
            split[int] - split ploly line to number of sections
            cut_dir['cw'|'ccw'] - cutting dir for closed polygons

        global (profile level):
            start - indicates the starting point of the cut

        gcode making steps:
        1. segments - all curves in dxf layer are convetred to lines, which in general, are unsorted. Those lines (2 3d points) are called segments later. A collection of segments form several layers make an unsorted profile.
        2. transformation - in order to sort segments and convert to a single polyline, following transformations are applied:
            *move to local CSYS
            *move to coresponding cross-section plane
        3. polyline - sorted cross section segments
        4. profile - polylines transformed into corss section planes p0, p1
        5. projection to axis - profiles projected to axes a0, a1
        6. cut model - group of prjected profiles, organize according to layer naming scheme
        7. generate gcode
    """
    conf = {'Z_span': 480, 'RTable_loc': 240}

    fname_dxf = args.input
    lname_dxf = args.layer
#find layers with the pattern and sort according to the convention
    seq_list, seq_dict = layers2seq(fname_dxf, lname_dxf)

    if seq_list and seq_dict:
    #CuttingSpace is a single cutting sequence
        ct = CuttingSpace(conf)
    #ModelSpace is a group of cut sequences to obtain the model
        md = ModelProfile()

        for j,  seq_name in enumerate(seq_list):
    #AxisProfile is an axis cutting sequence
            ap = AxisProfile()

            for i, axis_name in enumerate(seq_name):
    #create chain object
                s1 = chain(fname_dxf)

                for prof_name in axis_name:
    #add segments to chain. In general segments can be split between diferent layers
                    s1.AddSeg(prof_name)
                    # s1.MakeSplit()
    #apply transformations to segments like translation and rotation
                s1.ApplyTransformations()
    #sort segments to obtain a polyline
                s1.Seg2Poly()
                # s1.PlotChain()
    #add polylines to axis
                ap.Add2Axis(i, s1)
                # ap.Plot(mode='2D')
    #add axis to the model space
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
    parser.add_argument('-i', '--input', type=str, help='input filenames', required = True)
    parser.add_argument('-l', '--layer', type=str, help='input layers', required = True)
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
    # print(args)

    main(args)
