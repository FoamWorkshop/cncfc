#!/usr/bin/python
import numpy as np
import argparse
import matplotlib.pyplot as plt

def distance(p1,p2):
    # print(p1, p2)
    return np.linalg.norm(p2-p1)

def find_min_dist_point(arr,p):
    xy_dist=np.apply_along_axis(distance,1,arr,p)

    # print(xy_dist)
    return np.argmin(xy_dist)

parser = argparse.ArgumentParser(description='points to g code')
parser.add_argument('-i', '--input', type=str, required=True, help='input filename')
parser.add_argument('-zn', '--neutral_plane', type=float, default=0, help='neutral plane height for transition movements')
parser.add_argument('-zc', '--cut_plane', type=float, default=-1, help='cutting depth in 1 step')
parser.add_argument('-s', '--sort', type=bool, default=True, help='sort to minimize distances between points')

args = parser.parse_args()

f_in_name = args.input

zn = args.neutral_plane
zc = args.cut_plane
sort = args.sort

f_out_name = '{}.gcode'.format(f_in_name)
f_in=open(f_in_name,'r')
coords=f_in.readlines()
f_in.close()

xy_str = np.array([coord.split() for coord in coords])
xy=xy_str[:,:2].astype(float).round(1)
# print(xy)
pool = xy.copy()
# idx_pool=np.ones_like(xy)*None
# print(idx_pool)
p0=np.array([0,0])
xy_sorted=np.zeros_like(xy)
if sort:
    print('sort')
    # print(pool)
    for i, p in enumerate(xy):
        idx=find_min_dist_point(pool,p0)
        # print(idx)
        # print('pool size:',pool.shape)
        xy_sorted[i] = pool[idx]
        p0 = pool[idx]
        # print(p0)
        pool=np.delete(pool,idx,0)

    xy_sorted=np.unique(xy_sorted,axis=0)
    print(xy)

plt.figure(1)
ax1=plt.subplot(211)
ax1.plot(xy[:,0],xy[:,1],'o-')
ax1=plt.subplot(212)
ax1.plot(xy_sorted[:,0],xy_sorted[:,1],'o-')
plt.show()
if sort:
    xy=xy_sorted
header = 'G21 G90 G40 F800'
sufix = 'M84'

z_mov = [zn, zc]

f_out=open(f_out_name,'w')

f_out.write('{}\n'.format(header))

f_out.write('G1 Z {}\n'.format(z_mov[0]))
f_out.write('G1 X {0} Y {1}\n'.format(0,0))
# f_out.write('G1 Z {}\n'.format(z_mov[1]))

for i, coord in enumerate(xy):
    f_out.write('G1 Z {}\n'.format(z_mov[0]))
    f_out.write('G1 X {0} Y {1}\n'.format(coord[0],coord[1]))
    f_out.write('G1 Z {}\n'.format(z_mov[1]))

f_out.write('G1 Z {}\n'.format(z_mov[0]))
f_out.write('{}\n'.format(sufix))

f_out.close()
