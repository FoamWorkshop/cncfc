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
