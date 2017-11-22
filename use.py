#! /usr/bin/env python
from __future__ import absolute_import, division, print_function
import os
import argparse
import prog_path

overallPath = prog_path.pyPlanetPath
constitPath = prog_path.constituentPath
constituents = os.listdir(constitPath)

o = argparse.ArgumentParser(prefix_chars='-+')
for c in constituents:
    cpth = os.path.join(constitPath, c)
    if os.path.isdir(cpth):
        acon = '--' + c
        o.add_argument(acon, help='toggle ' + c, action='store_true')
args = o.parse_args()
argdict = vars(args)

for c in constituents:
    cpth = os.path.join(constitPath, c)
    if os.path.isdir(cpth):
        fil = os.listdir(cpth)
        if 'use.txt' in fil:
            fp = open(os.path.join(cpth, 'use.txt'), 'r')
            s = fp.readline().strip()
            fp.close()
            if argdict[c]:
                print('--- {:8s} Toggle to not use {}'.format(c, s))
                os.rename(os.path.join(cpth, 'use.txt'), os.path.join(cpth, 'nouse.txt'))
            else:
                print('+++ {:8s} Using {}'.format(c, s))
        elif 'nouse.txt' in fil:
            fp = open(os.path.join(cpth, 'nouse.txt'), 'r')
            s = fp.readline().strip()
            fp.close()
            if argdict[c]:
                print('+++ {:8s} Toggle to use {}'.format(c, s))
                os.rename(os.path.join(cpth, 'nouse.txt'), os.path.join(cpth, 'use.txt'))
            else:
                print('--- {:8s} Not using {}'.format(c, s))
        else:
            print('    {:8s} Use/nouse.txt not found'.format(c))
