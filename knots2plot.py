#!/usr/bin/python
'''the purpose of the program is to generate 3d model for FreeCAD from 2 given sections'''
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
col_list = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
arrow_length = 30
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

        line, = ax.plot(x1, y1, lw=2, marker=plt_marker,
                        color=col_list[i], label=input_file)
        lines.append(line)
        seg_len = ((x1[1] - x1[0])**2 + (y1[1] - y1[0])**2)**0.5
        len_rat = arrow_length / seg_len
        x2 = x1[0] + (x1[1] - x1[0]) * len_rat
        y2 = y1[0] + (y1[1] - y1[0]) * len_rat
        ax.annotate('({0:.2f},{1:.2f})'.format(x1[0], y1[0]),
                    xy=(x2, y2),
                    xytext=(x1[0], y1[0]),
                    arrowprops=dict(fc=col_list[i], ec=col_list[
                                    i], width=4, frac=0.1, shrink=0.1),
                    fontsize=6,
                    horizontalalignment='center',
                    verticalalignment='center',)
# head_width=0.05, head_length=0.1,
        if i == len(col_list) - 1:
            i = 0
        else:
            i += 1

        # arrows1.append(ax.plot(x1[0], y1[0],lw=1, marker='s',color='r',ms=10)[0])
        # arrows2.append(ax.plot(x1[1], y1[1],lw=1, marker='*',color='r',ms=10)[0])
        # for no, (x, y) in enumerate(zip(x1,y1)):
        #     # print('{0} ({1}, {2})'.format(i,x,y))
        # ax.text(x,y,'{0}'.format(no+1,x,y),fontsize=6) # ({1:.2f}, {2:.2f})

# legend display
leg = ax.legend(loc='upper left')
leg.get_frame().set_alpha(0.4)
# legend display

lined = dict()


for legline, origline in zip(leg.get_lines(), lines):
    legline.set_picker(5)  # 7 pts tolerance
    lined[legline] = origline


def onpick(event):
    # on the pick event, find the orig line corresponding to the
    # legend proxy line, and toggle the visibility
    legline = event.artist
    origline = lined[legline]
    vis = not origline.get_visible()
    origline.set_visible(vis)
    # Change the alpha on the line in the legend so we can see what lines
    # have been toggled
    if vis:
        legline.set_alpha(1.0)
    else:
        legline.set_alpha(0.2)
    fig.canvas.draw()

fig.canvas.mpl_connect('pick_event', onpick)

plt.axis('equal')
plt.grid(True)
plt.show()
