import numpy as np
from stl import mesh
import cncfclib as cf

# def main():
"""
This program converts stl files to paths for CNCFoamCutter.
slender models eg. fusselage modeled in OpenVSP.
"""

mesh = mesh.Mesh.from_file('fuselage.stl')
mesh.rotate([0,1,0],np.pi/2)

v_arr = np.round(np.vstack(mesh.vectors).astype(float), decimals=1)
#pos = [r, th, z] sectionwise
pos = cf.cartesian2cylyndrical(v_arr,17)
pos = cf.add_pedestal(pos)

profiles=np.zeros_like(pos)
for i in np.arange(np.shape(pos)[0]):
    profiles[i] = cf.cylyndrical2cartesian(pos[i])

    # cf.plot_loft_paths(profiles)
    # cf.plot_loft_paths(pos)

strokes = np.flipud(np.rot90(profiles))
# cf.plot_loft_paths(strokes)
# print(strokes)
ang_arr, r_arr, z_arr, v_arr = cf.transform(strokes, add_pedestal_bottom=True,add_pedestal_top=True)
cf.plot_surf(ang_arr,z_arr,r_arr)


res_dict = {'a_arr':ang_arr,
            'r_arr':r_arr,
            'z_arr':z_arr,
            'v_arr': v_arr}

np.save('res_dict.npy', res_dict)

print('saved')

# if __name__ == "__main__":
#     main()
