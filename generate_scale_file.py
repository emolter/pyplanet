from __future__ import print_function, absolute_import, division

import sys
import atmosphere


def generate(planet, constituent_list, value_list, calc_type='generate_calc', mode='mcmc', config='config.par', output_filename='Scratch/scale.dat', include_P=True):
    """
    Only need to include those constituents that change.
    Can include pressure (include_P) or not (just helpful for plotting - it gets ignored)
    Parameters:
    ------------
    planet:  planet name to find right atmosphere
    constituent_list:  list of constituent names that change, not case sensitive
    value_list:  xxx
    mode:  state_variable mode
    output_filename:  name where the data get written
    include_P:  boolean to add in pressure of not (convenience for plotting, etc)
    """
    sys.path.append(planet.capitalize())
    __import__(calc_type)
    generate_calc = sys.modules[calc_type]

    atm = atmosphere.Atmosphere(planet, mode=mode, config=config, log=None, plot=False)
    atm.run()
    atmospheric_pressure = atm.gas[atm.config.C['P']]
    atmospheric_temperature = atm.gas[atm.config.C['T']]
    atmospheric_value = {}
    for c in constituent_list:
        atmospheric_value[c] = atm.gas[atm.config.C[c.upper()]]

    with open(output_filename, 'w') as fp:
        cl = ' '.join(constituent_list)
        if include_P:
            s = '#p ' + cl + '\n'
        else:
            s = '#' + cl + '\n'
        fp.write(s)
        for i, pt in enumerate(zip(atmospheric_pressure, atmospheric_temperature)):
            s = ''
            if include_P:
                s = '{} '.format(pt[0])
            for c, v in zip(constituent_list, value_list):
                q = atmospheric_value[c][i]
                s += generate_calc.get_value(pt[0], pt[1], q, c, v) + ' '
            s = s.strip() + '\n'
            fp.write(s)
