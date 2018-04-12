from __future__ import absolute_import, division, print_function
import string
from scipy.interpolate import interp1d
import numpy as np
import prog_path
import atmosphere
import chemistry


def regrid(atm, regridType=None, Pmin=None, Pmax=None):
    """This puts atm and cloud on the same grid used later for calculations.  It has three different regimes:
        1 - interpolating within given points:  currently linear -- moving to: P,T,z are assumed to follow given adiabat
        2 - extrapolating outward:  project last slope out for everything
        3 - extrapolating inward:  uses a dry adiabat
       It has two regridTypes:
        1 - pressure grid points in file
            format: 'filename' (string)
                where
                    'filename' is a file with the pressures in bars (increasing pressure)
        2 - fixed number of steps in logP
            format:  'number'
                    where 'number' is number of steps within range (int)
    """

    verbose = False
    # set regridType/regrid
    if regridType is None:
        regridType = atm.config.regridType
    if (isinstance(regridType, str) and string.lower(regridType) == 'none') or regridType is None:
        print('No regridding.  Note that there is a risk that not everything is on the same grid...\n')
        return 0

    # set default Pmin/Pmax
    if Pmin is None or Pmin == 'auto' or Pmin == 0:
        Pmin = min(atm.gas[atm.config.C['P']])
    if Pmax is None or Pmin == 'auto' or Pmax == 0:
        Pmax = max(atm.gas[atm.config.C['P']])

    # set Pgrid
    if isinstance(regridType, str):
        try:
            regridType = int(regridType)
        except ValueError:
            Pgrid = _procF_(os.path.join(atm.config.path, regridType))
    if isinstance(regridType, int):
        regrid_numsteps = regridType
        Pgrid = np.logspace(np.log10(Pmin), np.log10(Pmax), regrid_numsteps)
        if verbose:
            print('Regridding on {} steps'.format(regrid_numsteps))
    Pmin = min(Pgrid)
    Pmax = max(Pgrid)
    # ## Copy over for new array size
    nAtm = len(Pgrid)
    nGas = atm.gas.shape[0]
    gas = np.zeros((nGas, nAtm))
    nCloud = atm.cloud.shape[0]
    cloud = np.zeros((nCloud, nAtm))

    # ## Interpolate gas onto the grid
    fillval = -999.9
    gas = interpolate('gas', gas, fillval, atm, Pgrid)

    # ## Extrapolate gas if needed
    if fillval in gas:
        atm.computeProp(False)
        gas = extrapolate(gas, fillval, atm)
    atm.gas = gas

    # ## Interpolate cloud onto the grid - fillval=0.0 extrapolates since we assume no clouds outside range and
    # ##      don't care about other stuff then either
    fillval = 0.0
    atm.cloud = interpolate('cloud', cloud, fillval, atm, Pgrid)

    # renormalize such that zDeep = 0.0 and reset DZ
    atm.renorm_z('gas')
    atm.renorm_z('cloud')

    print('Regrid:  {} levels'.format(nAtm))
    return 1


def _procF_(filename):
    try:
        reg = open(filename, 'r')
    except IOError:
        print("file '" + filename + "' not found - no regrid")
        return 0
    print('Regridding on file ' + filename)
    Pgrid = []
    for line in reg:
        Pgrid.append(float(line))
    reg.close()
    Pgrid = np.array(Pgrid)
    if not np.all(np.diff(Pgrid) > 0.0):
        print('Error in regrid:  P not increasing - flipping around and hoping for the best')
        Pgrid = np.fliplr(Pgrid)
    return Pgrid


def interpolate(gctype, gas_or_cloud, fillval, atm, Pgrid):
    # ## Interpolate gas/cloud onto the grid - currently both linear
    berr = False
    interpType = 'linear'

    if gctype == 'gas':
        ind = atm.config.C
        atm_gc = atm.gas
    else:
        ind = atm.config.Cl
        atm_gc = atm.cloud
    Pinput = atm_gc[ind['P']]

    gas_or_cloud[ind['P']] = Pgrid
    for yvar in ind:
        if yvar in ['P', 'DZ']:
            continue
        fv = interp1d(Pinput, atm_gc[ind[yvar]], kind=interpType, fill_value=fillval, bounds_error=berr)
        gas_or_cloud[ind[yvar]] = fv(Pgrid)

    return gas_or_cloud


def extrapolate(gas, fillval, atm):
    """First extrapolate in, then the rest of the fillvals get extrapolated out"""
    # extrapolate constituents in as fixed mixing ratios
    for yvar in atm.config.C:
        if yvar in ['Z', 'P', 'T', 'DZ']:
            continue
        val = atm.gas[atm.config.C[yvar]][-1]
        for i, V in enumerate(gas[atm.config.C[yvar]]):
            if V == fillval:
                gas[atm.config.C[yvar]][i] = val
    # extrapolate T and z as dry adiabat in hydrostatic equilibrium
    pDeep = atm.gas[atm.config.C['P']][-1]
    g = atm.layerProperty[atm.config.LP['g']][-1]
    r = atm.layerProperty[atm.config.LP['R']][-1]
    for i, p in enumerate(gas[atm.config.C['P']]):
        if p < pDeep:
            continue
        prev = gas[atm.config.C['P']][i - 1]
        dP = p - prev
        T = gas[atm.config.C['T']][i - 1]
        z = gas[atm.config.C['Z']][i - 1]
        amu = 0.0
        cp = 0.0
        for key in atm.chem:
            cp += atm.chem[key].cp * gas[atm.config.C[key]][i]
            amu += atm.chem[key].amu * gas[atm.config.C[key]][i]
        dT = T / (cp * p) * dP
        gas[atm.config.C['T']][i] = T + dT

        g = g + 2.0 * chemistry.R * T * np.log(p / prev) / (r * amu) / 1000.0
        H = chemistry.R * T / (amu * g) / 1000.0
        dz = H * dP / p
        r = r - dz
        gas[atm.config.C['Z']][i] = z - dz

    # extrapolate out along last slope (move return gas above if you think you want this)
    if fillval in gas:
        for yvar in atm.config.C:
            if yvar in ['P', 'DZ']:
                continue
            gas[atm.config.C[yvar]] = extrapolate_outward(gas[atm.config.C['P']], gas[atm.config.C[yvar]], fillval)

    return gas


def extrapolate_outward(x, y, fillval):
    for i in range(len(y)):
        if y[i] != fillval:
            break
    slope = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
    intercept = y[i] - slope * x[i]
    i = 0
    while y[i] == fillval:
        y[i] = slope * x[i] + intercept
        i += 1
    return y
