import sys
import os


def read_data(f_name, msg='False'):

    data = []

    if os.path.isfile(f_name):

        with open(f_name, 'r') as f:
            for line in f:
                tmp = line.split()
                x, y, z = 0, 0, 0
                if len(tmp) == 2:
                    x, y = tmp
                    z = 0
                if len(tmp) == 3:
                    x, y, z = tmp
                data.append([float(x), float(y), float(z)])
        if msg:
            print("{0:<24} -> {1} knots".format(f_name, len(data)))

    else:
        if msg:
            print('{0} not found'.format(f_name))

    return data


def write_data(f_name, data, msg='False'):
    i = 0

    bak_f_name = f_name+'.00.bak'

    while os.path.isfile(bak_f_name):
        bak_f_name = '{0}.{1:{fill}>2s}.bak'.format(f_name, str(i), fill='0')
        i += 1

    os.rename(f_name, bak_f_name)
    if msg:
        print("{0:<24} -> {1}".format(f_name, bak_f_name))

    with open(f_name, 'w') as f:
        for line in data:
            x, y, z = line
            f.write('{0:.6f} {1:.6f} {2:.6f}\n'.format(
                float(x), float(y), float(z)))

    if msg:
        print("{0:<24} <- {1} knots".format(f_name, len(data)))
