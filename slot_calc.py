import dxfgrabber
import numpy as np
import argparse
from cncfclib import *


def p_l_intersection(p0, vec_n, l0, l1):
    vec_l = np.subtract(l1, l0)
    param1 = np.subtract(p0, l0)
    d = (np.dot(param1, vec_n)) / (np.dot(vec_l, vec_n))
    vec_l = np.multiply(d, vec_l)
    return np.add(vec_l, l0)

p0 = [0, 0, 100]
vec_n = [0, 0, 1]

XY_p0 = [0, 130, 0]
XY_vec_n = [0, 1, 0]

UV_p0 = [0, -130, 0]
UV_vec_n = [0, 1, 0]

A0 = [97,       6.35,   -130]
B0=	[100.48,	-3.01,	-130]


A1 = [73.88,    10.85,   130]
B1=	[76.3,	     4.2,	130]


A2=p_l_intersection(p0, vec_n, A0, A1)
B2=p_l_intersection(p0, vec_n, B0, B1)

XY0=p_l_intersection(XY_p0, XY_vec_n, A0, B0)
UV0=p_l_intersection(UV_p0, UV_vec_n, B0, A0)

XY2=p_l_intersection(XY_p0, XY_vec_n, A2, B2)
UV2=p_l_intersection(UV_p0, UV_vec_n, A2, B2)

print("A2: {0}".format(A2))
print("B2: {0}".format(B2))
print("XY0: {0}".format(XY0))
print("UV0: {0}".format(UV0))
print("XY2: {0}".format(XY2))
print("UV2: {0}".format(UV2))
