
def generate(planet, constituent_list, value_list, output_filename='Scratch/scale.dat'):
    atm_pressure = planet.atm.gas[planet.atm.config.C['P']]
    atm_temperature = planet.atm.gas[planet.atm.config.C['T']]
    atm_value = {}
    chem = {}
    for c in constituent_list:
        atm_value[c] = planet.atm.gas[planet.atm.config.C[c.upper()]]
        chem[c] = planet.atm.chem[c.upper()]

    new_values = {}
    for c, v in zip(constituent_list, value_list):
        new_v = []
        for p, t, q in zip(atm_pressure, atm_temperature, atm_value[c]):
            s = get_value(p, t, q, c, v, chem)
            new_v.append(s)
        new_values[c] = np.asarray(new_v)


def get_value(p, t, q, c, v, chem):
    if c == 'H2S':
        Psat_gas = chem[c].Psat(t)
        if p < 43. and p * q * v > Psat_gas:  # Pressure greater than saturation pressure
            return 1.0
        elif p < 43. and p * q * v < Psat_gas:
            return v
        else:
            return 1.0
