#!/usr/bin/python
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='points to g code')
parser.add_argument('-i', '--input', type=str, required=True, help='input filename')
parser.add_argument('-zn', '--neutral_plane', type=float, default=0, help='neutral plane height for transition movements')
parser.add_argument('-zc', '--cut_plane', type=float, default=-1, help='cutting depth in 1 step')

args = parser.parse_args()

f_in_name = args.input

zn = args.neutral_plane
zc = args.cut_plane

f_out_name = '{}.gcode'.format(f_in_name)
f_in=open(f_in_name,'r')

coords=f_in.readlines()
f_in.close()

xy = np.array([coord.split() for coord in coords])

header = 'G21 G90 G40 F800'
sufix = 'M84'

z_mov = [zn, zc]

f_out=open(f_out_name,'w')

f_out.write('{}\n'.format(header))

# f_out.write('G1 Z {}\n'.format(z_mov[0]))
# f_out.write('G1 X {0} Y {1}\n'.format(0,0))
# f_out.write('G1 Z {}\n'.format(z_mov[1]))

for i, coord in enumerate(xy):
    f_out.write('G1 Z {}\n'.format(z_mov[0]))
    f_out.write('G1 X {0} Y {1}\n'.format(coord[0],coord[1]))
    f_out.write('G1 Z {}\n'.format(z_mov[1]))

f_out.write('G1 Z {}\n'.format(z_mov[0]))
f_out.write('{}\n'.format(sufix))

f_out.close()
