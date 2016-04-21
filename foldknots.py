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

FREECADPATH = '/usr/lib/freecad/lib/' # path to your FreeCAD.so or FreeCAD.dll file
import sys
sys.path.append(FREECADPATH)
import FreeCAD
#import PartGui
from FreeCAD import Base
import Part
import os
import argparse
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

def cross(a, b):
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]

    return c

def knots2face(knots, t, offset):
    knots_list = knots [:]
    if knots_list[0] != knots_list[-1]:
        knots_list.append(knots_list[0])

    u1=knots_list[1][0]-knots_list[0][0]
    u2=knots_list[1][1]-knots_list[0][1]
    u3=knots_list[1][2]-knots_list[0][2]
    v1=knots_list[-2][0]- knots_list[-1][0]
    v2=knots_list[-2][1]- knots_list[-1][1]
    v3=knots_list[-2][2]- knots_list[-1][2]

    if (u1*v1-u2*v2)>0:
        knots_list.reverse()

    u1=knots_list[1][0]-knots_list[0][0]
    u2=knots_list[1][1]-knots_list[0][1]
    u3=knots_list[1][2]-knots_list[0][2]
    v1=knots_list[-2][0]- knots_list[-1][0]
    v2=knots_list[-2][1]- knots_list[-1][1]
    v3=knots_list[-2][2]- knots_list[-1][2]

    points_list=[Base.Vector(var[0] , var[1] , var[2]) for var in knots_list]

    poly=Part.makePolygon(points_list)
    face=Part.Face(poly)
    if t==0:
        return face
    else:
        norm=cross([u1,u2,u3],[v1,v2,v3])
        norm_len=(norm[0]**2+norm[1]**2+norm[2]**2)**0.5
        print 'normal vector',[i/norm_len for i in norm]
        return face.extrude(Base.Vector(norm[0]/norm_len*-t,norm[1]/norm_len*-t,norm[2]/norm_len*-t))
        #return face.extrude(norm)

def interp_series(x,y,p):
    if x[0]>=p:
        return interp([x[0],y[0],0],[x[1],y[1],0],p)[1]
    if x[-1]<p:
        return interp([x[-2],y[-2],0],[x[-1],y[-1],0],p)[1]
    for i in range(len(x)-1):
        if x[i]<p<=x[i+1]:
            return interp([x[i],y[i],0],[x[i+1],y[i+1],0],p)[1]
def dist(x1,y1,x2,y2):
    return ((x2-x1)**2+(y2-y1)**2)**0.5

def interp_len(x,y): #converts x coord to a path distance from starting point s
    path_points=zip(x,y)
    path_segments_point=zip(path_points,path_points[1::])
    path_segments_lengt=[dist(s1[0],s1[1],s2[0],s2[1]) for s1, s2 in path_segments_point]

    path_dist = 0.0
    path_segments = [0]
    for seg_len in path_segments_lengt:
        path_dist+=seg_len
        path_segments.append(path_dist)

    return path_segments



#*********************************************************************DEFAULT PARAMETERS
dflt_o_f = 'all'  #decimal accuracy
dflt_dec_acc = 3  #decimal accuracy
dflt_n_arc = 10   #number of segments
dflt_l_arc = 1    #minimal segment length
dflt_path_dir = 1 #closed path collecting direction
sheet_thickness = 0
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

prof_f = 'ct_sec2_0.knt'
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
# [3, 0,  0],
# [5, 0, 0],
# [10,5, 0],
# [10,10, 0],
# [5, 15, 0],
# [0, 15, 0],
# [0, 0, 0]
# ]
# fold_knots=[
# [0, 0, 0],
# [2, 0, 0],
# [4, 0, 0],
# [6, 0, 0],
# [8, 5, 0],
# [10,20, 0],
# ]
print('{0}'.format('='*20))

prof_Xmin=min([var[0] for var in prof_knots])
prof_Xmax=max([var[0] for var in prof_knots])
prof_Ymin=min([var[1] for var in prof_knots])
prof_Ymax=max([var[1] for var in prof_knots])

print('profile:\n{4}Xmin: {0}\n{4}Xmax: {1}\n{4}Ymin: {2}\n{4}Ymax: {3}'.format(prof_Xmin, prof_Xmax, prof_Ymin, prof_Ymax,' '*8))
folds_list=[var for var in fold_knots if prof_Xmin < var[0] < prof_Xmax]
n_folds=len(folds_list)
print('folds: {0}'.format(n_folds))

pool=[]
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
print '\ngenerating pool knots: DONE'

#=============OK section
print '\npool knots'
Vertex_list=[]
for var in pool:
    print('{0:6.2f} {1:6.2f} {2:6.2f}'.format(var[0], var[1], var[2]))

faces_list=[]
folds_X = []

folds_X.append(prof_Xmin)
folds_X.extend([var[0] for var in folds_list])
folds_X.append(prof_Xmax)
print '\n\n'
for b1, b2 in zip(folds_X,folds_X[1::]):
    faces_list.append([var for var in pool if b1 <= var[0] <= b2])

#App.getDocument('testowy').addObject('App::DocumentObjectGroup','folded')
#for i, face in enumerate(faces_list):
    # print face
#    myPart=doc.addObject('Part::Feature','face_plan{0}'.format(i))
#    myPart.Shape = knots2face(face)
# print '\nprofile knots'}
# for var in prof_knots:
#     print var
#
print '\nfold knots'

fold_faces_list=[]
X=[var[0] for var in fold_knots]
Y=[var[1] for var in fold_knots]

print 'start folding'

for face in faces_list:
    fold_faces_list.append([[p1[0],p1[1],interp_series(X,Y,p1[0])] for p1 in face])

doc=FreeCAD.newDocument('testowy')#'testowy.fcstd')
myPart=doc.addObject('Part::Feature','fold_sheet')
myPart.Shape=Part.Compound([knots2face(face,sheet_thickness,0) for face in fold_faces_list])

#=============
unfold_faces_list=[]

X=[var[0] for var in fold_knots]
Y=[var[1] for var in fold_knots]
acc=0
print '\nunfolded knots'
uface=[]
fold_len_list=[]
# #     print var

unfold_faces_list=[]
fc_unfold_faces_list=[]
X=[var[0] for var in fold_knots]
Y=[var[1] for var in fold_knots]
x_len=interp_len(X,Y)
x0_len = interp_series(X,x_len,pool[0][0])
print 'start unfolding'

unfld_grp=doc.addObject('App::DocumentObjectGroup','unfold_sheet')
for face in faces_list:
   fc_unfold_faces_list.append([[interp_series(X,x_len,p1[0])-x0_len,p1[1],0] for p1 in face])

myPart=doc.addObject('Part::Feature','unfold_sheet')
myPart.Shape=Part.Compound([knots2face(face,sheet_thickness,0) for face in fc_unfold_faces_list])
unfld_grp.addObject(myPart)

#for face in faces_list:
unfold_faces_list.append([[interp_series(X,x_len,p1[0])-x0_len,p1[1],0] for p1 in pool])


for i, face in enumerate(unfold_faces_list):
    print 'face ',i
    for var in face:
        print var

doc.saveAs('testowy.fcstd')
