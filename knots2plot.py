#!/usr/bin/python
'''the purpose of the program is to visualise knot sets'''
# import os
# import sys
import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as npl
import xml.etree.ElementTree as ET
import argparse


from cncfclib import *
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import interp1d
from scipy.spatial import distance



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



lines = []
arrows_list = []
coords_list = []
coords_no_list = []
col_list = ['b', 'g', 'r', 'c', 'm', 'y', '0.75']
styl_list = ['-', '-', ':', ':', '--', '--', '-.']
al = 30
hl = al * 0.5
parser = argparse.ArgumentParser(description='plot knots')
parser.add_argument('-i',
                    '--input',
                    nargs='+',
                    type=str,
                    default=[],
                    help='input knots filenames')
parser.add_argument('-m',
                    '--marker',
                    type=str,
                    default='',
                    help='marker type')

args = parser.parse_args()

input_file_list = args.input
plt_marker = args.marker

if not input_file_list:
    input_file_list.extend([i for i in os.listdir('.') if i.endswith('.knt')])
    input_file_list.sort()

fig, ax = plt.subplots()
i = 0
for input_file in input_file_list:

    data = read_data(input_file, msg=True)
    if data:
        x1 = [var[0] for var in data]
        y1 = [var[1] for var in data]
# plot paths
        line, = ax.plot(x1, y1, lw=2, marker=plt_marker,
                        color=col_list[i], label=input_file, linestyle=styl_list[i])
        lines.append(line)

# plot direction arrows
        seg_len = ((x1[1] - x1[0])**2 + (y1[1] - y1[0])**2)**0.5
        len_rat = al / seg_len
        x2 = x1[0] + (x1[1] - x1[0]) * len_rat
        y2 = y1[0] + (y1[1] - y1[0]) * len_rat
        #({0:.2f},{1:.2f})
        arrow = ax.annotate(''.format(x1[0], y1[0]),
                            xy=(x2, y2),
                            xytext=(x1[0], y1[0]),
                            arrowprops=dict(fc=col_list[i], ec=col_list[
                                            i], width=4, headlength=hl, shrink=0.1),
                            fontsize=6,
                            horizontalalignment='center',
                            verticalalignment='center',)
        arrows_list.append(arrow)
# plot knot numbers and coordinates
        for no, var in enumerate(zip(x1, y1)):

            coord_no = ax.annotate('{0}'.format(no), xy=(var[0], var[1]), xytext=(var[0], var[1]),fontsize=6,horizontalalignment='center',
                verticalalignment='center',)

            # coord = ax.annotate('({0:.2f}, {1:.2f})'.format(
            #     var[0], var[1]), xy=(var[0], var[1]), xytext=(var[0], var[1]),fontsize=6,horizontalalignment='center',
            #     verticalalignment='center',)

# move to the next color
        if i == len(col_list) - 1:
            i = 0
        else:
            i += 1

# legend display
leg = ax.legend(loc='upper left')
leg.get_frame().set_alpha(0.4)
# legend display

lined = dict()
arrowed = dict()

for legline, origline in zip(leg.get_lines(), lines):
    legline.set_picker(5)  # 7 pts tolerance
    lined[legline] = origline

for legline, origarrow in zip(leg.get_lines(), arrows_list):
    arrowed[legline] = origarrow


def onpick(event):
    # on the pick event, find the orig line corresponding to the
    # legend proxy line, and toggle the visibility

    legline = event.artist

    origline = lined[legline]
    origarrow = arrowed[legline]
    vis = not origline.get_visible()
    origline.set_visible(vis)
    origarrow.set_visible(vis)

    # Change the alpha on the line in the legend so we can see what lines
    # have been toggled
    if vis:
        legline.set_alpha(1.0)
    else:
        legline.set_alpha(0.2)
    fig.canvas.draw()

fig.canvas.mpl_connect('pick_event', onpick)

plt.axis('equal')
plt.grid(b=True, which='major', color='k', linestyle='-', alpha=0.1)
plt.grid(b=True, which='minor', color='r', linestyle='-', alpha=0.1)
plt.minorticks_on()
#plt.grid(b=True, which='minor')#, linestyle='--')
plt.show()
