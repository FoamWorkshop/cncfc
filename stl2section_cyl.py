import numpy as np
from stl import mesh
import cncfclib as cf

mesh = mesh.Mesh.from_file('fuselage.stl')
mesh.rotate([0,1,0],np.pi/2)

v_arr = np.round(np.vstack(mesh.vectors).astype(float), decimals=1)
pos = cf.cartesian2cylyndrical(v_arr,17)

profiles=np.zeros_like(pos)
for i in np.arange(np.shape(pos)[0]):
    profiles[i] = cf.cylyndrical2cartesian(pos[i])

# cf.plot_loft_paths(profiles)

# print(profiles)

strokes = np.flipud(np.rot90(profiles))
# print(strokes)
# data=np.zeros_like(profiles)

# for i in np.arange(np.shape(profiles)[1]):
#     strokes[i]=profiles[:,i]
#     # ax.plot(spars[:,0],spars[:,1],spars[:,2],'x-')
#     ax.plot(data[:,i,0],data[:,i,1],data[:,i,2],'x-')
#
# print(profiles)
# b=np.flipud(np.rot90(a))

cf.interp_points_cyl(strokes,'sdsd','sadada')
