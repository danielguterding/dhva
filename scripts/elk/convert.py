import sys
from math import ceil
import numpy as np


def get_num_kpoints(lines):
    return [int(n) for n in lines[7].strip().split()]


def get_kvecs_and_remove_2Pi(lines):
    kvecs = []
    for l in lines[9:12]:
        kvecs.append([float(k) * 0.5 / np.pi for k in l.strip().split()])
    return kvecs


def write_header(fin, fout, lines, nx, ny, nz, kvecs):
    fout.writelines(lines[:7])
    fout.write(' {} {} {}\n'.format(nx - 1, ny - 1, nz - 1))
    fout.writelines(lines[8:9])
    for kx, ky, kz in kvecs:
        fout.write(
            '     {0: 1.8f}     {1: 1.8f}     {2: 1.8f}\n'.format(kx, ky, kz))
    fout.writelines(lines[12:13])


def get_energies_and_convert(lines, nx, ny, nz, startidx=13):
    energies = []
    ntot = nx * ny * nz
    for i, l in enumerate(lines[startidx:], startidx):
        energies += [float(e) for e in l.strip().split()]
        if len(energies) >= ntot:
            break
    energies = np.array(energies)
    energies.resize((nx - 1, ny - 1, nz - 1))
    energies *= 2  # Hartree to Rydberg, assumes Fermi Energy is zero
    energies = energies.flatten()
    return energies, i + 1


def write_energies(fout, energies, entries_per_line=6):
    num_lines = ceil(len(energies) / entries_per_line)
    for i in range(num_lines):
        dnidx = i * entries_per_line
        upidx = min(len(energies), dnidx + entries_per_line)
        elstr = ' ' + ' '.join(['{0: 1.6e}'.format(e)
                               for e in energies[dnidx:upidx]])
        fout.write(elstr + '\n')


def write_footer(fout, lines, endidx):
    fout.writelines(lines[endidx:])


if __name__ == '__main__':
    if not 3 == len(sys.argv):
        raise Exception('Please supply input and output file names.')
    infilename = sys.argv[1]
    outfilename = sys.argv[2]

    with open(infilename, 'rt') as fin:
        with open(outfilename, 'wt') as fout:
            lines = fin.readlines()
            nx, ny, nz = get_num_kpoints(lines)
            kvecs = get_kvecs_and_remove_2Pi(lines)
            write_header(fin, fout, lines, nx, ny, nz, kvecs)
            energies, endidx = get_energies_and_convert(lines, nx, ny, nz)
            write_energies(fout, energies)
            write_footer(fout, lines, endidx)
