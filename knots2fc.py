'''knots2fc:
-i1, i2, i3 ...
-solid i1, i2
-ruled i1, i2
-monge i1, i2
-t
-z

options:
  -i1, i2 ,i3... generate faces limited by input knots
  -solid: if only -i1 specified generate extruded solid with thickness z 
        if -i1 and -i2 specified generate lofted solid through 2 sections: 
            if -t then the distance between sections is t, otherwise, uses knots Z coordinates

  -monge: if only -i1 specified generate surface with i1 profile and  with thickness t and z distance between profiles
        -i1 XY projection -i2 XZ projection
TBD:
  -ruled: if only -i1 specified generate surface with i1 profile and  with thickness t and z distance between profiles
'''
