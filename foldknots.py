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

1. make list of profile sections
2. make list of fold knots fitting between selected sections

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

def cross_prod(u, v):
    return u[0] * v[1] - u[1] * v[0]

def sub_points(p1,p2):
    vect=[]
    p1=[x for x in p1[0]]
    p2=[x for x in p2[0]]
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

def knots2file(name,io_path,sorted_knots):
    f=open(name,'w')

    for var in io_path:
        coord=knot2coord(sorted_knots,var[0])
        f.write('{0:.2f} {1:.2f}\n'.format(coord[0],coord[1]))

    coord=knot2coord(sorted_knots, io_path[-1][1])
    f.write('{0:.2f} {1:.2f}\n'.format(coord[0],coord[1]))

    f.close()

def list_entities(dxf):
    dxf_summary=[shape.dxftype for shape in dxf.entities]
    print('{0:<10}: {1}'.format('LINES',dxf_summary.count('LINE')))
    print('{0:<10}: {1}'.format('ARCS',dxf_summary.count('ARC')))

def knots_dict(knots_list):
    return [[i,var] for i, var in enumerate(list(set(knots_list)))]

def elements_coords2knots(el_list, kt_list):
    el_kt_list=[]
    for el in el_list:
        for kt in kt_list:
            if kt[1]==el[0]: p1=kt[0]
            if kt[1]==el[1]: p2=kt[0]
        el_kt_list.append([p1,p2])
    return el_kt_list

def elements_knots2coords(el_list, kt_list):
    el_coord_list=[]
    for el in el_list:
        for kt in kt_list:
            if kt[0]==el[0]: p1=kt[1]
            if kt[0]==el[1]: p2=kt[1]
        el_coord_list.append([p1,p2])
    return el_coord_list

def knots_rank_list_summary(knots_rank):
    print('{0:<16}: {1}'.format('IO knots', [x[1] for x in knots_rank].count(1),\
                                        [x[0] for x in knots_rank if x[1]==1]))
    print('{0:<16}: {1}'.format('master knots',[x[1] for x in knots_rank].count(3),\
                                        [x[0] for x in knots_rank if x[1]==3]))
    print('{0:<16}: {1}'.format('chain knots', [x[1] for x in knots_rank].count(2),\
                                        [x[0] for x in knots_rank if x[1]==2]))

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

def find_path(crit,el_kt_list, sorted_knots, excl_knot):
    path=[]
    bar_length=50
    buff_length=len(el_kt_list)
    sect_length=int(2*buff_length/bar_length)
    counter=0
    knots_rank=knots_rank_list(el_kt_list, sorted_knots, excl_knot)
    val_max   =max(knots_rank, key=lambda tup: tup[1])[1]
#    print('number of elements: {0}'.format(buff_length))
#    print('0|{0}|100%'.format('-'*bar_length))
#    print ' |',
    while val_max>crit:
        knots_rank=knots_rank_list(el_kt_list, sorted_knots, excl_knot)
        curr_knot =knots_rank_find(knots_rank,1)
        val_max   =max(knots_rank, key=lambda tup: tup[1])[1]
        for element in el_kt_list:
            if curr_knot[0] in element:
                counter+=1
         #       print 'found element', element, 'val max', val_max
                if element[0]==curr_knot[0]:
                    path.append(element)
                if element[1]==curr_knot[0]:
                    path.append(element[::-1])

                el_kt_list.remove(element)
                break
            if counter>sect_length:
 #               print'+',
                counter=0
 #   print '+|DONE\n'
    if crit==1:
        path.append([path[-1][1],path[0][0]])
    return path

def find_l_el(read_dir,el_kt_list, sorted_knots, master_knot):
    #find all elements including master_knot and put into el_index list
    el_index=[i for i, element in enumerate(el_kt_list) if master_knot in element]

    seg1, seg2 = elements_knots2coords([el_kt_list[i] for i in el_index],sorted_knots)

    if cw_order(seg1, seg2)==read_dir: cur_ind=el_index[1]
    else: cur_ind=el_index[0]

    last_el=el_kt_list.pop(cur_ind)
    excl_knot=[x for x in last_el if x!=master_knot] #take the other knot

    return (last_el, excl_knot)

def cw_order(seg1,seg2):
    common_el=[x for x in list(set(seg1) & set(seg2))]
    u= sub_points(common_el,list(set(seg1)-set(common_el)))
    v= sub_points(common_el,list(set(seg2)-set(common_el)))
    if cross_prod(u, v)>0:
        return False
    else:
        return True

def print_list(data_list, common_text):
    for var in data_list:
        print('{0} {1}'.format(common_text, var))

def dxf_read(files, layer_name, dec_acc, n_arc, l_arc):
    tol=dec_acc
    knots_list=[]
    elements_list=[]
    line_count = 0
    arc_count = 0
#    list_entities(dxf)

    for shape in dxf.entities:
        if shape.layer==layer_name:
            if shape.dxftype=='LINE':
                line_count+=1
                p1=tuple(round(x,tol) for x in shape.start)
                p2=tuple(round(x,tol) for x in shape.end)
                knots_list.append(p1)
                knots_list.append(p2)
                elements_list.append([p1,p2])

            if shape.dxftype=='ARC':
                arc_count+=1
                ARC_knots_list=[]
                n=n_arc#number of segments
                min_len=l_arc
                O=shape.center
                R=shape.radius
                angl_1=shape.startangle*np.pi/180
                angl_2=shape.endangle*np.pi/180

                if angl_2>=angl_1:
                    angl_list=np.linspace(angl_1, angl_2, n)
                else:
                    angl_list=np.linspace(angl_1, angl_2+2*np.pi, n)

                arc_len=R*np.absolute(angl_2-angl_1)

                if arc_len/n<min_len:
                    n=max(int(arc_len/min_len),3)

                for angl in angl_list:
                    ARC_knots_list.append((round(O[0]+R*np.cos(angl),tol),round(O[1]+R*np.sin(angl),tol),O[2]))

                for i in range(n-1):
                    elements_list.append(ARC_knots_list[i:i+2])

                knots_list.extend(ARC_knots_list)

    return (knots_list, elements_list, [line_count, arc_count])

def save_knots(file_name,knots):
    f=open(file_name,'w')

    for var in knots:
        f.write('{0:.2f} {1:.2f}\n'.format(var[0], var[1]))

    f.close()

def read_knots(f_name):
    f=open(f_name,'r')
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
    f.close()

    print("input: {0} knots: {1}".format(f_name,len(data)))
    return data

def interp(s1, s2, p):
    s1x, s1y = s1[:-1]
    s2x, s2y = s2[:-1]
    return [p, (s2y - s1y) / (s2x - s1x) * (p - s1x) + s1y, 0]

# def interp(s1,s2,p):
#     x1, y1 = s1[:-1]
#     x2, y2 = s2[:-1]
#     return [p, (y2 - y1)/(x2 - x1)*(p-x1) + y1, 0]

def interp_series(x,y,p):
    if x[0]>=p:
        return interp([x[0],y[0],0],[x[1],y[1],0],p)[1]
    if x[-1]<p:
        return interp([x[-2],y[-2],0],[x[-1],y[-1],0],p)[1]
    for i in range(len(x)-1):
        if x[i]<p<=x[i+1]:
            return interp([x[i],y[i],0],[x[i+1],y[i+1],0],p)[1]



#*********************************************************************DEFAULT PARAMETERS
dflt_o_f = 'all'  #decimal accuracy
dflt_dec_acc = 3  #decimal accuracy
dflt_n_arc = 10   #number of segments
dflt_l_arc = 1    #minimal segment length
dflt_path_dir = 1 #closed path collecting direction
#*********************************************************************PROGRAM

# parser = argparse.ArgumentParser(description='test')
# parser.add_argument('-i1','--prof',tyoe=str,required=True,help='top view projection')
# parser.add_argument('-i2','--fold',type=str,required=True,help='side view projection')
# parser.add_argument('-o','--output',type=str,help='output file name')
#
# args = parser.parse_args()

# prof_f = args.prof
# fold_f = args.fold
# o_f = args.output

prof_f = 'i_sec2_0.knt'
fold_f = 'i_sec2_prof.knt'

prof_knots = read_knots(prof_f)
fold_knots = read_knots(fold_f)

o_list = []

# for s1, s2 in zip(prof_knots,prof_knots[1::]):
#     print s1, s2
#     fold_candidates=[ var for var in fold_knots if var[0] >= s1[0] and var[0] <= s2[0]]
#     print 'fold knots candidates',fold_candidates
#     # if fold_candidates:
#     #     for coord in fold_candidates:
#     #         o_list.append()
#     # fold_candidates = []
#
# print interp([1,1,0],[2,8,0],1.15)

# prof_knots=[
# [0, 0, 0],
# [3, 3,  0],
# [5, 0, 0],
# [10,5, 0]
# ]
# fold_knots=[
# [0, 0, 0],
# [2, 2, 0],
# [4, 2, 0],
# [6, 0, 0],
# [8, 0, 0],
# [10,3, 0],
# ]
print('{0}'.format('='*20))

print '\nprofile knots'
for var in prof_knots:
    print var

print '\nfold knots'
for var in fold_knots:
    print var
pool=[]
print '\n'

#=============OK section
for s1, s2 in zip(prof_knots,prof_knots[1::]):
    fit_candidates=[var[0] for var in fold_knots if s1[0]<var[0]<s2[0] or s1[0]>var[0]>s2[0]]
    print s1, s2
    pool.append(s1)
    if fit_candidates:
        if s2[0] >= s1[0]:
            for extra_point in fit_candidates:
                point = interp(s1,s2,extra_point)
                pool.append(point)
        else:
            for extra_point in fit_candidates[::-1]:
                point = interp(s2,s1,extra_point)
                pool.append(point)


pool.append(prof_knots[-1])

#=============OK section
print '\npool knots'
for var in pool:
    print('{0:6.2f} {1:6.2f} {2:6.2f}'.format(var[0], var[1], var[2]))

fold_sheet=[]
X=[var[0] for var in fold_knots]
Y=[var[1] for var in fold_knots]

for p1 in pool:
    fold_sheet.append([p1[0],p1[1],interp_series(X,Y,p1[0])])

print '\nfolded sheet'
for var in fold_sheet:
    print('{0:6.2f} {1:6.2f} {2:6.2f}'.format(var[0], var[1], var[2]))

#=============
unfold_sheet=[]

X=[var[0] for var in fold_knots]
Y=[var[1] for var in fold_knots]

for s1, s2 in zip(pool,pool[1::]):
    x2, y2 = s2[:-1]
    x1, y1 = s1[:-1]

    fold_y = interp_series(X,Y,x2)
    fold_len = (fold_y**2 + x2**2)**0.5
    unfold_sheet.append([x1 + fold_len, y2, 0])


print '\nunfolded knots'
for var in unfold_sheet:
    print('{0:6.2f} {1:6.2f} {2:6.2f}'.format(var[0], var[1], var[2]))






























# dir_path =os.getcwd()
# dxf_files=[i for i in os.listdir(dir_path) if i.endswith('.dxf')]
#
# if dxf_list != 'all':
#     print dxf_files
#     files_dxf=[i for i in dxf_files if i in dxf_list]
#
# if not files_dxf:
#     print 'dir does not include dxf files'
#
# #else, execute the program
# else:
#
#     print('SETTINGS:')
#     print('{0}{1:<30}: {2}'.format(' '*10,'decimal accuracy',dec_acc))
#     print('{0}{1:<30}: {2}'.format(' '*10,'arc segments count',n_arc))
#     print('{0}{1:<30}: {2}'.format(' '*10,'minimal arc segment length',l_arc))
#     print('{0}{1:<30}: {2}'.format(' '*10,'closed path collection dir',path_dir))
#     print('{0}{1:<30}: {2}'.format(' '*10,'files',files_dxf))
# #{{{print table header
#     print('\n\n{0:12}|{1:8}|{2:8}|{3:8}|{4:8}|{5:8}|{6:8}|{7:^20}'.format('file','layer','lines','arcs','1 -knt','2 -knt','3 -knt','status'))
#     print('{0}'.format('-'*80))
# #}}}print table header
#     for i, files_dxf_member in enumerate(files_dxf):
#
#         case_name=os.path.splitext(files_dxf_member)
#         dxf = dxfgrabber.readfile(files_dxf_member,{"assure_3d_coords": True})
#         dxf_layers = dxf.layers
#
#         for dxf_layers_member in dxf_layers:
#             layer_name = dxf_layers_member.name
#             knots_list, elements_list, shape_count = dxf_read(dxf, dxf_layers_member.name,dec_acc,n_arc, l_arc)
#             sorted_knots=knots_dict(knots_list)
#             el_kt_list  =elements_coords2knots(elements_list,sorted_knots)
#          #   print_list(el_kt_list, ' el kt list ')
#          #   print_list(sorted_knots,' sorted knots')
#             knots_rank  =knots_rank_list(el_kt_list, sorted_knots,None)
#             master_knot =knots_rank_find(knots_rank,3)
#             IO_knot     =knots_rank_find(knots_rank,1)
#
#             print('{0:12}|{1:8}|{2:8}|{3:8}|{4:8}|{5:8}|{6:8}|'.format(files_dxf_member,layer_name,\
#                                                                        shape_count[0],shape_count[1],\
#                                         [x[1] for x in knots_rank].count(1),\
#                                         [x[1] for x in knots_rank].count(2),\
#                                         [x[1] for x in knots_rank].count(3),' ')),
#
#
#             if len(IO_knot)!=1 or len(master_knot)!=1 or IO_knot[0]==None or master_knot[0]==None:
#                 print('{0:^20}|'.format('SKIPPED'))
#
#             else:
#
#                 io_path=find_path(2, el_kt_list, sorted_knots, None) #IO path
#                 print('{0:3}: {1:4d}|'.format('i/o',len(io_path))),
#
#                 last_el, excl_knot = find_l_el(path_dir, el_kt_list, sorted_knots, master_knot[0])
#
#                 ct_path=find_path(1, el_kt_list, sorted_knots, excl_knot[0]) #loop path
#                 print('{0:3}: {1:4d}|'.format('ct',len(ct_path)))
#
# #{{{
#                 i_file_name='{1}_{0}_{3}.{2}'.format(case_name[0],'i','knt',layer_name)
#                 knots2file(i_file_name, io_path, sorted_knots)
# #                print('i_path saved to: {0}'.format(i_file_name))
#
#                 o_file_name='{1}_{0}_{3}.{2}'.format(case_name[0],'o','knt',layer_name)
#                 knots2file(o_file_name, [ var[::-1] for var in io_path[::-1]] , sorted_knots)
# #                print('o_path saved to: {0}'.format(o_file_name))
#
#                 ct_file_name='{1}_{0}_{3}.{2}'.format(case_name[0],'ct','knt',layer_name)
#                 knots2file(ct_file_name, ct_path, sorted_knots)
# #                print('ct_path saved to: {0}'.format(ct_file_name))
# #}}}
