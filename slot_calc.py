from __future__ import division
import os
import sys
import numpy as np

def p_l_intersection(V0, n, P0, P1, coord='rel'):
    u = np.subtract(P1, P0)
    w = np.subtract(P0, V0)
    s = (-np.dot(n,w)) / (np.dot(n, u))
    su = np.multiply(s, u)
    
    Ps = np.add(w, su)

    if coord == 'abs':
        Ps = np.add(Ps , V0)
    
    return Ps

V0 = [0, 0, 100]
n = [0, 0, 1]

XY_V0 = [0, 130, 0]
XY_n = [0, 1, 0]

UV_V0 = [0, -130, 0]
UV_n = [0, 1, 0]

A0 =  [97,       6.35,   -130]
B0=	[100.48,	-3.01,	-130]


A1 = [73.88,    10.85,   130]
B1=	[76.3,	 4.2,	  130]


A2=p_l_intersection(V0, n, A0, A1, 'abs')
B2=p_l_intersection(V0, n, B0, B1, 'abs')

XY0=p_l_intersection(XY_V0, XY_n, A0, B0, 'abs')
UV0=p_l_intersection(UV_V0, UV_n, A0, B0, 'abs')

XY2=p_l_intersection(XY_V0, XY_n, A2, B2, 'abs')
UV2=p_l_intersection(UV_V0, UV_n, A2, B2, 'abs')

print("A2:  [{0[0]:6.2f}, {0[1]:6.2f}, {0[2]:6.2f}]".format(A2))
print("B2:  [{0[0]:6.2f}, {0[1]:6.2f}, {0[2]:6.2f}]".format(B2))
print("XY0: [{0[0]:6.2f}, {0[1]:6.2f}, {0[2]:6.2f}]".format(XY0))
print("UV0: [{0[0]:6.2f}, {0[1]:6.2f}, {0[2]:6.2f}]".format(UV0))
print("XY2: [{0[0]:6.2f}, {0[1]:6.2f}, {0[2]:6.2f}]".format(XY2))
print("UV2: [{0[0]:6.2f}, {0[1]:6.2f}, {0[2]:6.2f}]".format(UV2))

print ('results check:')
print ('XY0 UV0 vector length: {0:6.2f}'.format(np.linalg.norm(np.subtract(XY0,UV0))))
print ('XY UV planes distance: {0}'.format(260))
print ('angle: {0}'.format(np.rad2deg(np.arccos(260/np.linalg.norm(np.subtract(XY0,UV0)))))) 

print ('XY2 UV2 vector length: {0:6.2f}'.format(np.linalg.norm(np.subtract(XY2,UV2))))
print ('XY UV planes distance: {0}'.format(260))
print ('angle: {0}'.format(np.rad2deg(np.arccos(260/np.linalg.norm(np.subtract(XY2,UV2)))))) 
