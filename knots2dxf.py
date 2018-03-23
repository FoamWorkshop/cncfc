#!/usr/bin/python3
import ezdxf
import argparse
import re
def main(args):
    points=[]
    f_name = args.input
    with open(f_name,'r') as f:
        text = f.readlines()
        layer_name = text[0]
        for line in text[1:]:
            nodes = [ float(var) for var in re.findall('([\d\.-]+)',line)]
            print(nodes)
            points.append(nodes)


    dwg = ezdxf.new('AC1015')
    msp = dwg.modelspace()
    msp.add_lwpolyline(points)
    dwg.saveas(f_name + '.dxf')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='test')
    parser.add_argument('-i', '--input', help='airfoil dat profile')
    args = parser.parse_args()
    main(args)
