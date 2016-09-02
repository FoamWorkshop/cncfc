#!/usr/bin/python
'''the purpose of the program is to generate 3d model for FreeCAD from 2 given sections'''
FREECADPATH = '/usr/lib/freecad/lib/' # path to your FreeCAD.so or FreeCAD.dll file
import sys
sys.path.append(FREECADPATH)

import numpy as np
import numpy.linalg as npl
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt

import FreeCAD
import Part


from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import interp1d
from scipy.spatial import distance

#test

def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
lis    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    
    x_range = x_limits[1] - x_limits[0]; x_mean = np.mean(x_limits)
    y_range = y_limits[1] - y_limits[0]; y_mean = np.mean(y_limits)
    z_range = z_limits[1] - z_limits[0]; z_mean = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinityl
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_mean - plot_radius, x_mean + plot_radius])
    ax.set_ylim3d([y_mean - plot_radius, y_mean + plot_radius])
    ax.set_zlim3d([z_mean - plot_radius, z_mean + plot_radius])

def plot_path2d(data1, data2, data3, data4):
    x1=[x for [x, y, z] in data1]
    y1=[y for [x, y, z] in data1]
    z1=[z for [x, y, z] in data1]
    
    x2=[x for [x, y, z] in data2]
    y2=[y for [x, y, z] in data2]
    z2=[z for [x, y, z] in data2]

    mx1=[x for [x, y, z] in data3]
    my1=[y for [x, y, z] in data3]
    mz1=[z for [x, y, z] in data3]
    
    mx2=[x for [x, y, z] in data4]
    my2=[y for [x, y, z] in data4]
    mz2=[z for [x, y, z] in data4]

    pltxy=plt.plot(x1,y1, 'ro-',label='xy')
    pltuv=plt.plot(x2,y2, 'bs-',label='uv')
    mpltxy=plt.plot(mx1,my1, 'go-',label='cut_xy')
    mpltuv=plt.plot(mx2,my2, 'ks-',label='cut_uv')

    plt.legend()
    plt.axis('equal')
    plt.axis([min([min(x1),min(x2),min(mx1),min(mx2)]),\
              max([max(x1),max(x2),max(mx1),max(mx2)]),\
              min([min(y1),min(y2),min(my1),min(my2)]),\
              max([max(y1),max(y2),max(my1),max(my2)])]) 
    plt.grid(True)
    plt.interactive(True)
    plt.show(block=False)
#   plt.show()
 
    plt.hold(True)

def plot_path3d(data1, data2, data3, data4):
    x1=[x for [x, y, z] in data1]
    y1=[y for [x, y, z] in data1]
    z1=[z for [x, y, z] in data1]
    
    x2=[x for [x, y, z] in data2]
    y2=[y for [x, y, z] in data2]
    z2=[z for [x, y, z] in data2]

    mx1=[x for [x, y, z] in data3]
    my1=[y for [x, y, z] in data3]
    mz1=[z for [x, y, z] in data3]
    
    mx2=[x for [x, y, z] in data4]
    my2=[y for [x, y, z] in data4]
    mz2=[z for [x, y, z] in data4]
  
    fig=plt.figure()
    ax=fig.gca(projection='3d')
    pltxy=ax.plot(x1, y1, z1,'ro-',label='xy')
    pltuv=ax.plot(x2, y2, z2,'bs-',label='uv')
    mashpltxy=ax.plot(mx1, my1, mz1,label='mxy')
    mashpltuv=ax.plot(mx2, my2, mz2,label='muv')

    ax.legend()
    set_axes_equal(ax)
    plt.show() 
    plt.hold(True)

def p_l_intersection(p0,vec_n,l0,l1):
    vec_l=np.subtract(l1,l0)
    param1=np.subtract(p0,l0)
    d=(np.dot(param1,vec_n))/(np.dot(vec_l,vec_n))
    vec_l=np.multiply(d,vec_l)
    return np.add(vec_l,l0)

def p_l_intersection_series(p0,vec_n,data1,data2):
    if len(data1)==len(data2):
        print "data ok"
        tmp=[]
        for i in range(len(data1)):
            l0=data1[i]
            l1=data2[i]
            print l0, l1
            tmp.append(p_l_intersection(p0,vec_n,l0,l1))
        return tmp
    else:
        return [0,0,0]
   
class model:
    def __init__(self, sect1, sect2):
        self.sect1=sect1
        self.sect2=sect2
    
class path_points:
    '''simple'''
    def __init__(self, data):
        self.data=data

    def rotate(self, rc_x, rc_y, angl_deg):
        angl_rad=np.deg2rad(angl_deg)
        sn=np.sin(angl_rad)
        cs=np.cos(angl_rad)
        self.data=[(  cs*(line[0]-rc_x)- sn*(line[1]-rc_y)+rc_x,\
                      sn*(line[0]-rc_x)+ cs*(line[1]-rc_y)+rc_y,\
                      line[2]) for line in self.data]
 
    def translate(self, tr_x, tr_y, tr_z):
        self.data=[(line[0]+tr_x, line[1]+tr_y, line[2]+tr_z) for line in self.data]
 
    def scale(self, sc_x, sc_y):
        self.data=[(line[0]*sc_x, line[1]*sc_y, line[2]) for line in self.data]
        
def gcodexyuv(dataxy, datauv):

    if len(dataxy)==len(datauv):
        print "data ok"
        tmp=[]
        for i in range(len(dataxy)):
            tmp.append('g1 x{0:6.2f} y{1:6.2f} u{2:6.2f} v{3:6.2f}'.format(dataxy[i][0], dataxy[i][1], datauv[i][0], datauv[i][1]))
            
        fgcode=open('test.ngc','w')
        for line in tmp:
            fgcode.write(line)
            fgcode.write('\n')
        fgcode.close()
        return tmp

    else:
        print "nie mozna wygenerowac g codu. rozne dlugosci sciezek."
        return 0
    
def read_data(f_name):
    f=open(f_name,'r')
    data=[]
    
    for line in f:
        tmp=line.split()
        x,y,z=0,0,0
        if len(tmp)==2:
            x,y=tmp
            z=0
        if len(tmp)==3:
            x,y,z=tmp
        data.append([float(x),float(y),float(z)])
    f.close()
    print("\nsection: {0}    nodes: {1}".format(f_name,len(data)))
    return data


def PolyArea(x,y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

def silentremove(filename):
    if os.path.exists(filename):
        os.remove(filename)
def points2face(points):
    fc_profile1_elems=[]

    for p1, p2 in zip(points, points[1:]):
        fc_profile1_elems.append(Part.makeLine(p1, p2))

    fc_profile1_wire=Part.Wire(fc_profile1_elems)
    return Part.Face(fc_profile1_wire)

def balance_points(data1, data2):
    '''function balances 2 data sets, data2 to data1'''
    dist_1=[]
    dist_2=[]
    buff=[]
    x=[]
    y=[]
    z=[]
    data1_u=[]
    data1_N=len(data1)
    data2_N=len(data2)
    cos_uv=[]
    #  length_set=[]
    dist_1.append(0)
    dist_1.extend([distance.euclidean(u,v) for u, v in zip(data1, data1[1:])])
    length_1=np.cumsum(dist_1)
    
    dist_2.append(0)
    dist_2.extend([distance.euclidean(u,v) for u, v in zip(data2, data2[1:])])
    length_2=np.cumsum(dist_2)
    #print data1[6][0]
    data1_A=0
    
    for i in range(data1_N-1):
         data1_A+=0.5*( data1[i][0]*data1[i+1][1]-data1[i][1]*data1[i+1][0])

    data1_C1=0
    data1_C2=0
    for i in range(data1_N-1):
         data1_C1+=1/(6*data1_A)*(data1[i][0]+data1[i+1][0])*(data1[i][0]*data1[i+1][1]-data1[i][1]*data1[i+1][0])
         data1_C2+=1/(6*data1_A)*(data1[i][1]+data1[i+1][1])*(data1[i][0]*data1[i+1][1]-data1[i][1]*data1[i+1][0])
    
    print("total length data1: {0},\narea: {1},\nC1: {2}\nC2: {3}".format(length_1[-1], data1_A,data1_C1,data1_C2))       
    
    data1_u=[(var[0]-data1_C1,var[1]-data1_C2) for var in data1]
    data2_u=[(var[0]-data1_C1,var[1]-data1_C2) for var in data2]

    print("\n ====data1 vectors===")

    for var in data1_u:
        print var

    print("\n ====data2 vectors===")

    for var in data2_u:
        print var

    
    # data1_x=[ var[0] for var in data1]
    # data1_y=[ var[1] for var in data1]
    # data1_z=[ var[2] for var in data1]
   
   # print data1_x
    
    
#    f21=length_1[-1] / length_2[-1]
    
 #   print("total length data1: {0},\ntotal length data2: {1},\nconversion factor: {2}\n".format(length_1[-1], length_2[-1],f21))
    
    # length_set=list(set(length_1)|set(length_2))
    # length_set.sort()
    # for i in length_set:
    #     print i
    # print 'length 1'
    # for i in length_1:
    #     print i

    # print 'length 2'  
    # for i in length_2:
    #     print i

    # print 'length 2 *f21'  
    # for i in length_2:
    #     print i*f21

        
    # buff=[(np.interp(dl,length_1,data1_x),\
    #        np.interp(dl,length_1,data1_y),\
    #        np.interp(dl,length_1 ,data1_z)) for dl in length_set]

    
    # return buff

dist_z=2
span=500
dihedral=10
test_data1=[]

#dataxy=read_data('naca2416.dat')
#datauv=read_data('naca0012.dat')

datauv=read_data('naca0012.dat')
dataxy=read_data('naca2416.dat')
#test_data1=balance_points(dataxy,datauv)
#test_data2=balance_points(datauv,dataxy)

# print("\nproject datauv to dataxy")
# for var in test_data1:
#     print("{0}".format(var))

# print("\nproject dataxy to datauv")    
# for var in test_data2:
#     print("{0}".format(var))



    
pathxy=path_points(dataxy)
pathuv=path_points(datauv)

pathxy.scale(100, 100)
pathxy.rotate(0, 0, 0)
pathxy.translate(0, 0, 0.5*(dist_z+span))

pathuv.scale(60, 60)
pathuv.rotate(0, 0, 5)
pnewathuv.translate(0, 0, 0.5*(dist_z-span))

cut_obj=model(pathxy.data,pathuv.data)

mashpathxy=path_points(p_l_intersection_series([0,0,dist_z],[0,0,1],pathxy.data,pathuv.data))

mashpathuv=path_points(p_l_intersection_series([0,0,0],[0,0,1],pathuv.data,pathxy.data))

#for var in collection:
    
mashpathxy.scale(1, 1)
mashpathxy.rotate(0, 0, 0)
mashpathxy.translate(0, 0, dist_z)

mashpathuv.scale(0.7, 0.7)
mashpathuv.rotate(0, 0, 5)
mashpathuv.translate(10, 0, 0)

#plot_path2d(pathxy.data, pathuv.data, mashpathxy.data, mashpathuv.data)        
#plot_path3d(pathxy.data, pathuv.data, mashpathxy.data, mashpathuv.data)        

# # #gcodexyuv(pathxy.data, pathuv.data)


doc=FreeCAD.newDocument()#'testowy.fcstd')

myPart=doc.addObject('Part::Feature','wing')

#cube = Part.makeBox(2,2,2)
#p1=Part.Point(1,2,3)
#p2=Part.Point(123,311,312)
points1=pathxy.data
points2=pathuv.data

face1=points2face(points1)
face2=points2face(points2)
# fc_profile1_elems=[]

# for p1, p2 in zip(points1,points1[1:]):
#     fc_profile1_elems.append(Part.makeLine(p1, p2))

# fc_profile1_wire=Part.Wire(seg) fc_profile1_face=Part.Face(fc_profile1_wire)
poly1=Part.makePolygon(points1)
poly2=Part.makePolygon(points2)

myPart.Shape = Part.makeLoft([poly1,poly2])

#fc_profile1_face

#myPart.Shape = cube
#vo.show
#doc.recompute()
silentremove('test_data1')
doc.saveAs('testowy.fcstd')
#FreeCAD.closeDocument('test_data1.fcstd')

