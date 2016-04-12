#!/usr/bin/python

'''dxf2path.py, 6.04.2016 author: Adam Narozniak
dxf2path ptogram is a part of the CNCFCgcode generator. It automaticaly extracts cutting path from a dxf file. 
The cutting path is split into:
1. IO_path - in/out path begining with single knot
2. ct_path - closed loop, begin and end in master knot position
the output of the program is a set of files with an ordered list of knots'''

import os
import getopt
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
    print('number of elements: {0}'.format(buff_length))
    print('0|{0}|100%'.format('-'*bar_length))
    print ' |',
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
                print'+',
                counter=0
    print '+|DONE\n'           
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
  
def dxf_read(files,layer_name):
    tol=3
    knots_list=[]
    elements_list=[]

    list_entities(dxf)

    for shape in dxf.entities:
        if shape.layer==layer_name:
            if shape.dxftype=='LINE':
                p1=tuple(round(x,tol) for x in shape.start)
                p2=tuple(round(x,tol) for x in shape.end)
                knots_list.append(p1)
                knots_list.append(p2)
                elements_list.append([p1,p2])

            if shape.dxftype=='ARC':
                ARC_knots_list=[]
                n=20 #number of segments
                min_len=0.5
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

    return (knots_list, elements_list)

def save_knots(file_name,knots):
    f=open(file_name,'w')

    for var in knots:
        f.write('{0:.2f} {1:.2f}\n'.format(var[0], var[1]))

    f.close()


#****************************************************************************
#*****************************program**************************************** 
#****************************************************************************
dir_path=os.getcwd()
dxf_files=os.listdir(dir_path)

arg_len=len(sys.argv)
if arg_len == 1:
    print 'go through all dxfs in the folder\n'
    #find dxf files in the current directory
    files_dxf=[i for i in dxf_files if i.endswith('.dxf')]

else:
    cmdargs = sys.argv
    print cmdargs[1:]
    files_dxf=[i for i in dxf_files if i in cmdargs[1:]]

    #if list is empty display msg and exit 
if not files_dxf:

    print 'dir does not include dxf files' 

#else, execute the program

else:
    #list found files:
   # print '\n'
    print_list(files_dxf,'found: ')
  #  print '\n'
    for i, files_dxf_member in enumerate(files_dxf):

        case_name=os.path.splitext(files_dxf_member)

      
        dxf = dxfgrabber.readfile(files_dxf_member,{"assure_3d_coords": True})
        dxf_layers = dxf.layers

        for dxf_layers_member in dxf_layers:            
            layer_name = dxf_layers_member.name
            print('\n{0}'.format('#'*60))
            print('file {0}: {1}'.format(i, files_dxf_member))
            print('LAYER: {0}'.format(layer_name))
            print('{0}\n'.format('#'*60))                                 
            knots_list, elements_list = dxf_read(dxf, dxf_layers_member.name)
            sorted_knots=knots_dict(knots_list)
            el_kt_list  =elements_coords2knots(elements_list,sorted_knots)
            knots_rank  =knots_rank_list(el_kt_list, sorted_knots,None)
            master_knot =knots_rank_find(knots_rank,3)  
            IO_knot     =knots_rank_find(knots_rank,1)  

            knots_rank_list_summary(knots_rank)

            if len(IO_knot)!=1 or len(master_knot)!=1 or IO_knot[0]==None or master_knot[0]==None:
                print('{0}{1}'.format('-'*24,'SKIPPED'))

            else:
                print '\nIO path:'
                io_path=find_path(2, el_kt_list, sorted_knots, None) #IO path
            #    print 'IO path: DONE'
                last_el, excl_knot = find_l_el(1, el_kt_list, sorted_knots, master_knot[0])
                print '\nct path:'
                ct_path=find_path(1, el_kt_list, sorted_knots, excl_knot[0]) #loop path
            #    print 'ct path: DONE'
            #   paths_summary(io_path,ct_path)

                io_file_name='{1}_{0}_{3}.{2}'.format(case_name[0],'io','knt',layer_name)
                knots2file(io_file_name, io_path, sorted_knots)
                print('IO_path saved to: {0}'.format(io_file_name))

                ct_file_name='{1}_{0}_{3}.{2}'.format(case_name[0],'ct','knt',layer_name)
                knots2file(ct_file_name, ct_path, sorted_knots)
                print('ct_path saved to: {0}'.format(ct_file_name))
