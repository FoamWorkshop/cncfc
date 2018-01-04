import numpy as np
from stl import mesh
import cncfclib as cf

mesh = mesh.Mesh.from_file('fuselage.stl')
mesh.rotate([0,1,0],np.pi/2)

v_arr = np.round(np.vstack(mesh.vectors).astype(float), decimals=1)
pos = cf.cartesian2cylyndrical(v_arr,17)

spars=np.zeros_like(pos)
for i in np.arange(np.shape(pos)[0]):
    spars[i] = cf.cylyndrical2cartesian(pos[i])

cf.plot_loft_paths(spars)
