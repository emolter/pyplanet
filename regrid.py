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
    atm.computeProp(False)  # We need to do this to interpolate/extrapolate along the adiabat

    # ## Interpolate gas onto the grid
    fillval = -999.9
    gas = interpolate('gas', gas, fillval, atm, Pgrid)

    # ## Extrapolate gas if needed
    if fillval in gas:
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
    use_lapse_hydro = False  # This will use adiabat and hydrostatic for T and Z

    if gctype == 'gas':
        ind = atm.config.C
        atm_gc = atm.gas
        Pinput = atm_gc[ind['P']]
        if use_lapse_hydro:
            fv = interp1d(Pinput, atm.layerProperty[atm.config.LP['LAPSEP']], kind=interpType, fill_value=0.0, bounds_error=False)
            S = {}
            S['T'] = fv(Pgrid)  # Lapse rate
            fv = interp1d(Pinput, atm.layerProperty[atm.config.LP['H']], kind=interpType, fill_value=0.0, bounds_error=False)
            S['Z'] = fv(Pgrid) / Pgrid  # H/P
            dP = []
            xvar = {}
            xvar['Z'] = []
            xvar['T'] = []
            for p in Pgrid:
                for i, Ptrial in enumerate(Pinput):
                    if Ptrial > p:
                        break
                dP.append(Ptrial - p)
                xvar['Z'].append(atm.gas[atm.config.C['Z']][i])
                xvar['T'].append(atm.gas[atm.config.C['T']][i])
            dP = np.array(dP)
            testing = {}
    else:
        ind = atm.config.Cl
        atm_gc = atm.cloud
        Pinput = atm_gc[ind['P']]

    gas_or_cloud[ind['P']] = Pgrid
    for yvar in ind:
        if yvar in ['P', 'DZ']:
            continue
        if use_lapse_hydro and gctype == 'gas' and yvar in ['Z', 'T']:
            testing[yvar] = xvar[yvar] + dP * S[yvar]
            #gas_or_cloud[ind[yvar]] = xvar[yvar] + dP * S[yvar]
            #continue
        fv = interp1d(Pinput, atm_gc[ind[yvar]], kind=interpType, fill_value=fillval, bounds_error=berr)
        gas_or_cloud[ind[yvar]] = fv(Pgrid)

    if gctype == 'gas' and use_lapse_hydro:
        print("use_lapse_hydro does not work well and linear works really well")
        atmosphere.plt.figure("NEWZ")
        atmosphere.plt.plot(Pinput, atm_gc[ind['Z']], 'r.')
        atmosphere.plt.plot(Pgrid, testing['Z'], 'b+')
        atmosphere.plt.plot(Pgrid, gas_or_cloud[ind['Z']], 'gx')
        atmosphere.plt.figure("NEWT")
        atmosphere.plt.plot(Pinput, atm_gc[ind['T']], 'r.')
        atmosphere.plt.plot(Pgrid, testing['T'], 'b+')
        atmosphere.plt.plot(Pgrid, gas_or_cloud[ind['T']], 'gx')

    return gas_or_cloud


def extrapolate(gas, fillval, atm):
    """First extrapolate in, then the rest of the fillvals get extrapolated out"""
    Pmin = min(atm.gas[atm.config.C['P']])
    Pmax = max(atm.gas[atm.config.C['P']])
    # extrapolate constituents in as fixed mixing ratios
    for yvar in atm.config.C:
        if yvar in ['Z', 'P', 'T', 'DZ']:
            continue
        val = atm.gas[atm.config.C[yvar]][-1]
        for i, P in enumerate(gas[atm.config.C['P']]):
            if P > Pmax:
                gas[atm.config.C[yvar]][i] = val
    # extrapolate T and z as dry adiabat in hydrostatic equilibrium
    gDeep = atm.layerProperty[atm.config.LP['g']][-1]
    pDeep = atm.layerProperty[atm.config.LP['P']][-1]
    rDeep = atm.layerProperty[atm.config.LP['R']][-1]
    for i, P in enumerate(gas[atm.config.C['P']]):
        if P > Pmax:
            dP = P - gas[atm.config.C['P']][i - 1]
            P = gas[atm.config.C['P']][i - 1]
            T = gas[atm.config.C['T']][i - 1]
            amu = 0.0
            cp = 0.0
            for key in atm.chem:
                cp += atm.chem[key].specific_heat * gas[atm.config.C[key]][i]
                amu += atm.chem[key].amu * gas[atm.config.C[key]][i]
            cp *= chemistry.R  # since the catalogued values are Cp/R
            dT = (chemistry.R * T) / (cp * P) * dP
            gas[atm.config.C['T']][i] = T + dT

            g = gDeep + 2.0 * chemistry.R * T * np.log(P / pDeep) / (rDeep * amu) / 1000.0
            H = chemistry.R * T / (amu * g) / 1000.0
            dz = H * dP / P
            gas[atm.config.C['Z']][i] = gas[atm.config.C['Z']][i - 1] - dz
    return gas

    # extrapolate out along last slope (move return gas above if you think you want this)
    if fillval in gas:
        for yvar in atm.config.C:
            if yvar in ['P', 'DZ']:
                continue
            gas[atm.config.C[yvar]] = extrapolateOut(gas[atm.config.C['P']], gas[atm.config.C[yvar]], fillval)

    return gas


def extrapolateOut(x, y, fillval):
    if y[0] == fillval:
        slope = (y[b[0] + 1] - y[b[0]]) / (x[b[0] + 1] - x[b[0]])
        intercept = y[b[0]] - slope * x[b[0]]
        i = 0
        while y[i] == fillval:
            y[i] = slope * x[i] + intercept
            i += 1
    if y[-1] == fillval:
        slope = (y[b[-1]] - y[b[-1] - 1]) / (x[b[-1]] - x[b[-1] - 1])
        intercept = y[b[-1]] - slope * x[b[-1]]
        i = -1
        while y[i] == fillval:
            y[i] = slope * x[i] + intercept
            i -= 1
    return y
