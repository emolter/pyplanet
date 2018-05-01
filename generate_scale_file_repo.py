from __future__ import print_function, absolute_import, division

import atmosphere
import chemistry


def generate(planet, constituent_list, value_list, mode='mcmc', output_filename='Scratch/scale.dat', include_P=True):
    """
    Only need to include those constituents that change.
    Can include pressure (include_P) or not (just helpful for plotting - it gets ignored)
    Parameters:
    ------------
    planet:  planet name to find right atmosphere
    constituent_list:  list of constituent names that change, not case sensitive
    output_filename:  name where the data get written
    include_P:  boolean to add in pressure of not (convenience for plotting, etc)
    """
    atm = atmosphere.Atmosphere(planet, mode=mode, config='planet', log=None, plot=False)
    atm.run()
    atmospheric_pressure = atm.gas[atm.config.C['P']]
    atmospheric_temperature = atm.gas[atm.config.C['T']]

    with open(output_filename, 'w') as fp:
        cl = ' '.join(constituent_list)
        if include_P:
            s = '#p ' + cl + '\n'
        else:
            s = '#' + cl + '\n'
        fp.write(s)
        for p, t, q in zip(atmospheric_pressure, atmospheric_temperature, atmospheric_value):
            s = ''
            if include_P:
                s = '{} '.format(p)
            for c, v in zip(constituent_list, value_list):
                s += get_value(p, t, q, c, v) + ' '
            s = s.strip() + '\n'
            fp.write(s)


def get_value(p, t, q, c, v):
    """
    Put all the smart stuff here.  Can it be parameterized by a few good MCMC-able things?
    This is a sample from Josh...
    """

    gas = chemistry.ConstituentProperties(c)
    Psat_gas = gas.Psat(t)

    if c == 'H2S':
        if p < 43. and p * q * v > Psat_gas:  # Pressure greater than saturation pressure
            return str(1.0)
        elif p < 43. and p * q * v < Psat_gas:
            return str(v)
        else:
            return str(1.0)
