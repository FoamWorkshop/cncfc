import numpy as np
f_in=open('holes.txt','r')


coords=f_in.readlines()
f_in.close()

# for coord in coords:
#     i, j = coord.split()
#     print(float(x), float(y))

xy = np.array([coord.split() for coord in coords])
# .astype(float)

header = 'G21 G90 G40 F800'
sufix = 'G28'
# z_depth = -2
z_mov = [0, -2]

f_out=open('holes.gcode','w')
f_out.write('{}\n'.format(header))
for i, coord in enumerate(xy):
    f_out.write('G1 Z {}\n'.format(z_mov[0]))
    f_out.write('G1 X {0} Y {1}\n'.format(coord[0],coord[1]))
    f_out.write('G1 Z {}\n'.format(z_mov[1]))

f_out.write('G1 Z {}\n'.format(z_mov[0]))
f_out.write('{}\n'.format(sufix))

f_out.close()
