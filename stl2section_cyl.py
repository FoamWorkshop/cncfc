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
# cf.plot_loft_paths(pos)

strokes = np.flipud(np.rot90(profiles))
ang_arr, r_arr, z_arr, v_arr = cf.transform(strokes, add_pedestal_bottom=True,add_pedestal_top=True)

print('ang', ang_arr),
print('r', r_arr)
print('z', z_arr)
print('v', v_arr)

np.savetxt('r_arr',r_arr)
np.savetxt('a_arr',ang_arr)
np.savetxt('z_arr',z_arr)
# np.savetxt('v_arr',v_arr)
# add_pedestal()

# buff_2v.tofile('ndarray_test.dat',sep=";")
# np.save('ndarray_test.dat', buff_2v)
# k=np.load('ndarray_test.dat.npy', mmap_mode='r')
# print(k)
print('saved')
