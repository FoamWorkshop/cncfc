#!/usr/bin/python
import xml.etree.ElementTree as ET
import copy
import numpy as np
import argparse

# parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
#                                 description='''\
# --------------------------
# KiCAD svg PCB cleaner
# --------------------------
#
# 1. Removes white circles
# 2. Changes 0 length paths to circles
# 3. Changes circles size to 2mm
# ''',
#                                 usage='%(prog)s [-i dxf_filename] [-l lofted_shape_prefix]')


# scale=0.00254

scale=1
# /25.4/96


parser = argparse.ArgumentParser(description='clean svg')
parser.add_argument('-i', '--input', type=str, required=True, help='input filename')
parser.add_argument('-d', '--circle_diameter', type=float, default=2, help='input filename')

args = parser.parse_args()

f_in_name = args.input
d = args.circle_diameter

f_out_name = 'mod_{}'.format(f_in_name)
f_out_holes_name = 'holes_{}'.format(f_in_name)

ET.register_namespace('', 'http://purl.org/dc/elements/1.1/')
ET.register_namespace('', 'http://creativecommons.org/ns#')
ET.register_namespace('', 'http://www.w3.org/1999/02/22-rdf-syntax-ns#')
ET.register_namespace('', 'http://www.w3.org/2000/svg')

tree = ET.parse(f_in_name)

root = tree.getroot()
white_circle_list=[]
white_circle_cent=[]
path_0_list=[]

height = float(root.attrib['height'].split('mm')[0])
print(height)

for child in root:
    if 'circle' in child.tag and '#ffffff' in child.attrib['style']:
        white_circle_list.append(child)
        white_circle_cent.append([float(child.attrib['cx']),float(child.attrib['cy']),child.attrib['id']])

for child in root:
    if 'path' in child.tag and 'v 0' in child.attrib['d']:
        coords=child.attrib['d'].split()[1]
        x, y = coords.split(',')
        path_0_list.append([child,float(x),float(y)])
        print(float(x), float(y))

print('path 0: {}'.format(len(path_0_list)))
print('white circle: {}'.format(len(white_circle_list)))

for elm in white_circle_list:
    root.remove(elm)

r=d/2/scale

for i, (elm, x, y) in enumerate(path_0_list):
    root.remove(elm)
    circle_id = i + 9000
    new_circle = ET.Element('circle')
    new_circle.set('cx',str(x))
    new_circle.set('cy',str(y))
    new_circle.set('r',str(r))
    new_circle.set('id',str(circle_id))
    new_circle.set('style','fill:#0;fill-opacity:1;stroke:#840000;stroke-width:0;stroke-linecap:round;stroke-linejoin:round;stroke-opacity:1')
    root.append(new_circle)
    new_circle=[]

for child in root:
    if 'circle' in child.tag:
        child.set('r',str(r))
        child.set('style','fill:#0;fill-opacity:1;stroke:#840000;stroke-width:0;stroke-linecap:round;stroke-linejoin:round;stroke-opacity:1')

tree.write(f_out_name)



f=open(f_out_holes_name,'w')
for hole in white_circle_cent:
    f.write('{} {} {}\n'.format(hole[0]*scale, height-hole[1]*scale, hole[2]))
f.close()
