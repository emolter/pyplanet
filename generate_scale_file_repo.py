from __future__ import print_function, absolute_import, division

import atmosphere


def generate(planet, constituent_list, output_filename='Scratch/scale.dat', include_P=True):
    atm = atmosphere.Atmosphere(planet, config='planet', log=None, plot=False)
    atm.run()
    atmospheric_pressure = atm.gas[atm.config.C['P']]
    with open(output_filename, 'w') as fp:
        cl = ' '.join(constituent_list)
        if include_P:
            s = '#p ' + cl + '\n'
        else:
            s = '#' + cl + '\n'
        fp.write(s)
        for p in atmospheric_pressure:
            s = ''
            if include_P:
                s = '{} '.format(p)
            for c in constituent_list:
                s += get_value(p, c) + ' '
            s = s.strip() + '\n'
            fp.write(s)


def get_value(p, c):
    if c == 'nh3':
        v = str(0.5)
    else:
        v = str(0.9)
    return v
