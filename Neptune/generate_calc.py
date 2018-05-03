

def get_value(planet, p, t, q, c, v):
    """
    Put all the smart stuff here.  Can it be parameterized by a few good MCMC-able things?
    """

    gas = planet.atm.chem[c.upper()]
    Psat_gas = gas.Psat(t)

    if c == 'H2S':
        if p < 43. and p * q * v > Psat_gas: # Pressure greater than saturation pressure 
            return 1.0
        elif p < 43. and p * q * v < Psat_gas:
            return v
        else:
            return 1.0