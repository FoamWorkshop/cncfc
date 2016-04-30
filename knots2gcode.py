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

def cross_prod(u, v):
    return u[0] * v[1] - u[1] * v[0]

def sub_points(p1,p2):
    vect=[]
    p1=[x for x in p1[0]]
    p2=[x for x in p2[0]]
 #   print p3, p4
    print len(p1)
    print len(p2)
    if len(p1)==len(p2):
        for i, n in enumerate(p2):
            vect.append(n-p1[i])
        return vect
    return len(p1)*[None]

def knots_rank_find(knots_rank,rank):
    knots=[x[0] for x in knots_rank if x[1]==rank]
    if len(knots)>0: return knots
    else: return [None]

def knots_rank_list(el_kt_list,sorted_knots,skip_knot):
    knots_rank=[]
    for var in sorted_knots:
        if var[0]==skip_knot:
            knots_rank.append([var[0],None])
        else:
            knots_rank.append([var[0],[x for a in el_kt_list for x in a].count(var[0])])
    return knots_rank

def knot2coord(sorted_knots,knot):
    for var in sorted_knots:
        if var[0]==knot:
            return var[1]
    return None

def knots2gcode(ct_path,name='gcode',header='True'):
    with open(name,'w') as f:

        if header:
            f.write("G21 (Units in millimeters)\n")
            f.write("G90 (Absolute programming)\n")
            f.write("G40 (Cancel radius comp.)\n")
            f.write("G49 (Cancel length comp.)\n")

    ##############################################
        f.write("(-----CT PATHS-----{0}-----)\n".format(name))
        for var in ct_path:
            f.write('G1 X{0:.1f} Y{1:.1f} U{0:.1f} V{1:.1f}\n'.format(var[0], var[1]))
    ##############################################
        if header:
            f.write("M2 (Program end)")

    #f.close()

def read_data(f_name):
    with open(f_name) as f:
        data=[]

        for line in f:
            tmp=line.split()
            x,y,z=0,0,0
    #        print tmp
            if len(tmp)==2:
                z=0
                x, y = tmp
            if len(tmp)==3:
                x,y,z=tmp

            data.append([float(x),float(y),float(z)])


    print("input: {0} knots: {1}".format(f_name,len(data)))
    return data

# layer_name='cut'
# knots_list, elements_list = dxf_read(file_name+'.dxf',layer_name)
# sorted_knots=knots_dict(knots_list)
# el_kt_list  =elements_coords2knots(elements_list,sorted_knots)
# knots_rank  =knots_rank_list(el_kt_list, sorted_knots,None)
# master_knot =knots_rank_find(knots_rank,3)
# IO_knot     =knots_rank_find(knots_rank,1)
#
# knots_rank_list_summary(knots_rank)


#*********************************************************************DEFAULT PARAMETERS
dflt_dxf_list = 'all'  #decimal accuracy
dflt_dec_acc = 3  #decimal accuracy
dflt_n_arc = 10   #number of segments
dflt_l_arc = 1    #minimal segment length
dflt_path_dir = 1 #closed path collecting direction
#*********************************************************************PROGRAM
knt_data=[]
parser = argparse.ArgumentParser(description='test')
parser.add_argument('-i','--input',nargs='*',default=dflt_dxf_list,help='input filenames')

args = parser.parse_args()

knt_list = args.input

for knt_set in knt_list:
    print 'processing:', knt_set

    knt_data = read_data(knt_set)

    if knt_data:
        output_name = knt_set+'.ngc'
        knots2gcode(knt_data,output_name,False)
