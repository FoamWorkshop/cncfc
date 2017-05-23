import ezdxf

points=[]
f_name = 'naca4416.dat'
f=open(f_name,'r')

for line in f:
    line = line.strip(' ')
    line = line.strip('\n')
    line = line.split(',')

    if type(line)== str:
        line = line.split()
    # print()
    # print()
    points.append((float(line[0]), float(line[1])))
f.close()
print(points)
dwg = ezdxf.new('AC1015')
msp = dwg.modelspace()

# points = [(0, 0), (3, 0), (6, 3), (6, 6)]
# points = [(0, 0), (3, 0), (6, 3), (6, 6)]
msp.add_lwpolyline(points)
dwg.saveas(f_name + '.dxf')
