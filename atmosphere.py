from __future__ import absolute_import, division, print_function
import math
import string
import matplotlib.pyplot as plt
import numpy as np
import utils
import sys
import os
import os.path

# ##local imports
import prog_path
import config as pcfg
import raypath as ray
import properties
import regrid

planetDictionary = {'Jupiter': 0, 'Saturn': 1, 'Uranus': 2, 'Neptune': 3}


class Atmosphere:
    def __init__(self, planet, config='config.par', path=None, log=None, batch=False, verbose=False, plot=True):
        """reads/computes atmospheres.  This should return:
               self.gas
               self.cloud
               self.layerProperty
            on the appropriate grid
            Note that the research is in the input files and modifying the tweak modules
            All of the default config parameters are hard-coded here:  see __init__, setConfig, showConfig."""

        planet = string.capitalize(planet)
        self.planet = planet
        self.verbose = verbose
        self.plot = plot
        self.logFile = utils.setupLogFile(log)
        self.batch = batch

        print('\n---Atmosphere of {}---'.format(planet))
        if type(config) == str:
            config = pcfg.planetConfig(self.planet, configFile=config, log=log)
        self.config = config

        # ##Create function dictionaries
        self.gasGen = {}
        self.gasGen['read'] = self.readGas
        self.gasGen['compute'] = self.computeGas
        self.cloudGen = {}
        self.cloudGen['read'] = self.readCloud
        self.cloudGen['compute'] = self.computeCloud
        self.propGen = {}
        self.propGen['read'] = self.readProp
        self.propGen['compute'] = self.computeProp

        print('Planet ' + self.planet)
        if self.config.gasType == 'read':  # this assumes that cloudType is then also 'read'
            utils.log(self.logFile, '\tReading from: ' + self.config.filename, True)
            utils.log(self.logFile, '\tAtmosphere file:  ' + self.config.gasFile, True)
            utils.log(self.logFile, '\tCloud file:  ' + self.config.cloudFile, True)
        if verbose:
            print(self.config.show())

    def run(self, Pmin=None, Pmax=None, regridType=None, gasType=None, cloudType=None, otherType=None, tweak=True):
        """This is the standard pipeline"""
        # ##Set run defaults
        if Pmin is None:
            Pmin = self.config.pmin
        if Pmax is None:
            Pmax = self.config.pmax
        if regridType is None:
            regridType = self.config.regridType
        if gasType is None:
            gasType = self.config.gasType
        if cloudType is None:
            cloudType = self.config.cloudType
        if otherType is None:
            otherType = self.config.otherType
        self.nAtm = 0

        # ## Generate gas profile (gasType is 'read' or 'compute')
        if gasType not in self.gasGen.keys():
            print('Error:  No such gasType: ', gasType)
            return 0
        else:
            self.gasGen[gasType](verbose=self.verbose)

        if not self.batch:
            # ## Generate cloud profile (cloudType is 'read' or 'compute')
            if cloudType not in self.cloudGen.keys():
                print('Error:  No such cloudType: ', cloudType)
                return 0
            else:
                self.cloudGen[cloudType](verbose=self.verbose)

            if tweak:  # This loads and calls the module 'tweakFile'
                self.tweakAtm()

            # ## Compute other parameters that are needed
            if self.config.otherType not in self.propGen.keys():
                print('Error:  no such otherTpe: ', otherType)
                return 0
            else:
                self.propGen[otherType](verbose=self.verbose)

            # ## Put onto common grid
            regridded = regrid.regrid(self, regridType=regridType, Pmin=Pmin, Pmax=Pmax)
            self.nAtm = len(self.gas[0])

            angularDiameter = 2.0 * math.atan(self.layerProperty[self.config.LP['R']][0] / self.config.distance)
            if self.verbose:
                print('angular radius = {} arcsec'.format(utils.r2asec(angularDiameter / 2.0)))

            # ## Plot data
            if self.plot:
                self.plotTP()
                self.plotGas()
                self.plotCloud()
                self.plotProp()

        return self.nAtm

    def getval(self, val=None, vtype='gas'):
        """Returns one of the constituent or cloud profiles"""
        if val is None:
            print("Usage:  getVal('v',['gas'/'cloud'/'other'])")
            print("    'gas' is default")
            print('These are the gas values:')
            print(self.config.C)
            print('These are the cloud values:')
            print(self.config.Cl)
            print('These are the layerProperty (other) values:')
            print(self.config.LP)
            return None
        v = string.upper(val)
        vt = string.lower(vtype)
        rv = 0
        found = False

        if vt == 'gas':
            if v in self.config.C.keys():
                rv = self.gas[self.config.C[v]]
                print('Found ' + val + ' in gas')
                found = True
        elif vt == 'cloud':
            if v in self.config.Cl.keys():
                rv = self.cloud[self.config.Cl[v]]
                print('Found ', val, ' in cloud')
                found = True
        elif vt == 'other':
            if v in self.config.LP.keys():
                rv = self.layerProperty[self.config.LP[v]]
                print('Found ', val, ' in layerProperty')
                found = True
        print(vt, '  ', val, end='')
        if found:
            print(':  found')
        else:
            print(':  not found')
        return rv

    def plotTP(self, plot='auto'):
        """Plot the T-P profile"""
        if plot == 'auto':
            plt.figure(planetDictionary[self.planet])
        plt.title(self.planet + ':  T-P profile')
        plt.loglog(self.gas[self.config.C['T']], self.gas[self.config.C['P']])
        v = list(plt.axis())
        v[2] = 100.0 * math.ceil(self.gas[self.config.C['P']][-1] / 100.0)
        v[3] = 1.0E-7 * math.ceil(self.gas[self.config.C['P']][0] / 1E-7)
        plt.axis(v)
        plt.ylabel('P [bars]')
        plt.xlabel('T [K]')

    def plotCloud(self, dontPlot=['Z', 'P', 'T', 'DZ'], plot='auto'):
        """Plots the clouds"""
        if plot == 'auto':
            plt.figure(planetDictionary[self.planet] + 11)
        plt.title(self.planet + ': clouds')
        for cloud in self.config.Cl:
            present, cl = self.__isPresent__(self.cloud[self.config.Cl[cloud]])
            if cloud in dontPlot or not present:
                continue
            plt.loglog(cl, self.cloud[self.config.Cl['P']], label=cloud)
        v = list(plt.axis())
        if v[0] < 1E-10:
            v[0] = 1E-10
        v[2] = 100.0 * math.ceil(self.gas[self.config.C['P']][-1] / 100.0)
        v[3] = 1.0E-7 * math.ceil(self.gas[self.config.C['P']][0] / 1E-7)
        plt.axis(v)
        plt.ylabel('P [bars]')
        plt.xlabel(r'Density [g/cm$^3$]')
        plt.legend()

    def plotGas(self, dontPlot=['Z', 'P', 'T', 'DZ'], plot='auto'):
        """Plots the constituents"""
        if plot == 'auto':
            plt.figure(planetDictionary[self.planet] + 10)
        plt.title(self.planet + ': gas')
        for gas in self.config.C:
            present, g = self.__isPresent__(self.gas[self.config.C[gas]])
            if gas in dontPlot or not present:
                continue
            plt.loglog(g, self.gas[self.config.C['P']], label=gas)
        v = list(plt.axis())
        if v[0] < 1E-10:
            v[0] = 1E-10
        v[2] = 100.0 * math.ceil(self.gas[self.config.C['P']][-1] / 100.0)
        v[3] = 1.0E-7 * math.ceil(self.gas[self.config.C['P']][0] / 1E-7)
        plt.axis(v)
        plt.ylabel('P [bars]')
        plt.xlabel('Fractional Abundance')
        plt.legend()

    def plotProp(self, dontPlot=['Z', 'P', 'T'], plot='auto'):
        if plot == 'auto':
            plt.figure(planetDictionary[self.planet] + 12)
        plt.title(self.planet + ': other')
        for other in self.config.LP:
            present, g = self.__isPresent__(self.layerProperty[self.config.LP[other]])
            if other in dontPlot or not present:
                continue
            plt.loglog(g, self.gas[self.config.C['P']], label=other)
        v = list(plt.axis())
        if v[0] < 1E-10:
            v[0] = 1E-10
        v[2] = 100.0 * math.ceil(self.gas[self.config.C['P']][-1] / 100.0)
        v[3] = 1.0E-7 * math.ceil(self.gas[self.config.C['P']][0] / 1E-7)
        plt.axis(v)
        plt.ylabel('P [bars]')
        plt.xlabel('Property')
        plt.legend()

    def readGas(self, gasFile=None, Cdict=None, verbose=False):
        """Reads gas profile file as self.gas"""

        if gasFile is None:
            gasFile = self.config.gasFile
        if Cdict is None:
            Cdict = self.config.C
        if gasFile == 'batch':
            self.batch = True
            return 0
        else:
            self.batch = False
        gasFile = os.path.join(self.config.path, gasFile)

        print('Reading constituents from {}'.format(gasFile))
        self.gas = []
        print('\tUsing atmsopheric component:  ', end='')
        for k in Cdict:
            print(k + '  ', end='')
            self.gas.append([])
        self.nConstituent = len(Cdict)
        try:
            fp = open(gasFile, "r")
        except IOError:
            print(gasFile + ' was not found - returning no gas profile\n\n')
            raise IOError
        print(' ')
        for line in fp:
            if line[0] in utils.commentChars or len(line) < 4:
                continue
            data = line.split()
            try:
                cval = [float(x) for x in data]
            except ValueError:
                continue
            if len(cval) > self.nConstituent:
                print("{} has wrong number of constituents".format(line))
                continue
            for n in range(self.nConstituent):  # Initialize all of the constituents to 0.0
                if n < len(cval):
                    self.gas[n].append(cval[n])
                else:
                    self.gas[n].append(0.0)
        self.nGas = len(self.gas[0])

        # ##Redo the constituent dictionary for the self.gas index positions
        nid, sk = utils.invertDictionary(Cdict)
        for i, k in enumerate(sk):
            Cdict[nid[k]] = i
        self.config.C = Cdict

        print('\tRead ' + str(self.nGas) + ' lines')
        fp.close()
        self.gas = np.array(self.gas)
        # ## Check that P is monotonically increasing
        monotonic = np.all(np.diff(self.gas[self.config.C['P']]) > 0.0)
        if not monotonic:
            self.gas = np.fliplr(self.gas)
            monotonic = np.all(np.diff(self.gas[self.config.C['P']]) > 0.0)
        if not monotonic:
            print("Error in ", gasFile, ".  Pressure not monotonically increasing")

        # ## Renormalize so that deepest z is 0
        zDeep = self.gas[self.config.C['Z']][-1]
        if abs(zDeep) > 1E-6:
            for i in range(len(self.gas[self.config.C['Z']])):
                self.gas[self.config.C['Z']][i] -= zDeep
            print('Deepest levels set from %.1f to 0'.format(zDeep))

        # ## Set dz
        dz = np.abs(np.diff(self.gas[self.config.C['Z']])) * 1.0E5  # convert from km to cm (so no unit checking!!!)
        self.gas[self.config.C['DZ']] = np.append(np.array([0.0]), dz)
        return self.nGas

    def writeGas(self, outputFile='gas.dat'):
        fp = open(outputFile, 'w')
        gdat = np.transpose(self.gas)
        fp.write(str(len(gdat)) + '\n')
        for data in gdat:
            s = ''
            for i in range(11):
                s += (str(data[i]) + '\t')
            s += '\n'
            fp.write(s)
        fp.close()

    def readCloud(self, cloudFile=None, Cldict=None, verbose=False):
        """Reads in cloud data if we have it..."""

        if cloudFile is None:
            cloudFile = self.config.cloudFile
        if Cldict is None:
            Cldict = self.config.Cl
        if cloudFile == 'batch':
            self.batch = True
            return 0
        else:
            self.batch = False
        cloudFile = os.path.join(self.config.path, cloudFile)

        print('Reading clouds from {}'.format(cloudFile))
        self.cloud = []
        print('\tUsing cloud components:  ', end='')
        for k in Cldict:
            print(k, '  ', end='')
            self.cloud.append([])
        self.nParticulate = len(Cldict.keys())
        try:
            fp = open(cloudFile, "r")
        except IOError:
            print(cloudFile, ' was not found - returning no clouds\n\n')
            raise IOError
        print(' ')
        for line in fp:
            if line[0] in utils.commentChars or len(line) < 4:
                continue
            data = line.split()
            try:
                cval = [float(x) for x in data]
            except ValueError:
                continue
            if len(cval) > self.nParticulate:
                print("{} has wrong number of particulates".format(line))
                continue
            for n in range(self.nParticulate):  # Initialize all of the particulates to 0.0
                if n < len(cval):
                    self.cloud[n].append(cval[n])
                else:
                    self.cloud[n].append(0.0)
        self.nCloud = len(self.cloud[0])
        # ##Redo the particulate dictionary for the self.cloud index positions
        nid, sk = utils.invertDictionary(Cldict)
        for i, k in enumerate(sk):
            Cldict[nid[k]] = i
        self.config.Cl = Cldict

        print('\tRead ', str(self.nCloud), ' lines')
        fp.close()
        self.cloud = np.array(self.cloud)
        # ## Check that P is monotonically increasing
        monotonic = np.all(np.diff(self.cloud[self.config.Cl['P']]) > 0.0)
        if not monotonic:
            self.cloud = np.fliplr(self.cloud)
            monotonic = np.all(np.diff(self.cloud[self.config.Cl['P']]) > 0.0)
        if not monotonic:
            print("Error in ", cloudFile, ".  Pressure not monotonically increasing")
        self.cloud[self.config.Cl['DZ']] = np.abs(np.append(np.diff(self.cloud[self.config.Cl['Z']]), 0.0))
        return self.nCloud

    def readProp(self, otherFile=None, verbose=False):
        """Reads in other property data if we have it..."""
        if otherFile is None:
            otherFile = self.config.otherFile
        otherFile = os.path.join(self.config.path, otherFile)
        print("readProp currently doesn't read in any other file...")
        print("it will probably not be needed, since it will get computed in computeProp")
        return False

    def tweakAtm(self):
        """Tweaks the atmosphere data..."""
        if self.batch:
            print('Defer tweaking')
            return 0
        nAtm = len(self.gas[self.config.C['P']])
        if nAtm - self.nGas != 0:
            print('Error in number of layers - check it out (', str(nAtm), '/', str(self.nGas), ')')
            print('Returned from tweakAtm')
            return 0
        # Import tweakFile
        sys.path.append(self.config.path)
        try:
            __import__(self.config.tweakFile)
            tweakModule = sys.modules[self.config.tweakFile]
        except SyntaxError:
            utils.log(self.logFile, "Syntax Error:  check " + self.config.tweakFile, True)
            return 0

        # Run module then log
        self.tweakComment, self.gas, self.cloud = tweakModule.modify(self.gas, self.cloud, self.config.C, self.config.Cl)
        print('---tweakComment')
        print(self.tweakComment)
        print('---')
        utils.log(self.logFile, self.tweakComment, False)
        _tf = os.path.join(self.config.path, self.config.tweakFile + '.py')
        _tp = open(_tf, 'r')
        dt = _tp.read()
        utils.log(self.logFile, '======================' + _tf + '=====================', False)
        utils.log(self.logFile, dt, False)
        utils.log(self.logFile, '====================================================================', False)
        _tp.close()

        return nAtm

    def computeProp(self, verbose=False):
        """This module computes derived atmospheric properties (makes self.layerProperty)"""
        if self.batch:
            print('Defer computing properties')
            return 0
        nAtm = len(self.gas[self.config.C['P']])
        self.layerProperty = []
        for op in self.config.LP:
            self.layerProperty.append([])
        zOffset = 0.0
        iOffset = 0
        psep = 1.0E6
        for i, zv in enumerate(self.gas[self.config.C['Z']]):     # find the nearest z value at p_ref
            P = self.gas[self.config.C['P']][i]
            if abs(P - self.config.p_ref) < psep:
                psep = abs(P - self.config.p_ref)
                iOffset = i
        zOffset = self.gas[self.config.C['Z']][iOffset]
        z_at_p_ref = self.config.Req
        if verbose:
            print("z,P offset:  ", zOffset, self.gas[self.config.C['P']][iOffset])

        for i, zv in enumerate(self.gas[self.config.C['Z']]):
            T = self.gas[self.config.C['T']][i]
            P = self.gas[self.config.C['P']][i]
            self.layerProperty[self.config.LP['P']].append(P)
            self.layerProperty[self.config.LP['Z']].append(zv)
            rr = z_at_p_ref + zv - zOffset
            self.layerProperty[self.config.LP['R']].append(rr)   # note that this is the "actual" z along equator  referenced to planet center (aka radius)
            # ##set mean amu
            amulyr = 0.0
            for key in properties.amu:
                if key in self.config.C:
                    amulyr += properties.amu[key] * self.gas[self.config.C[key]][i]
            self.layerProperty[self.config.LP['AMU']].append(amulyr)
            # ##set GM pre-calc (normalized further down) and get lapse rate
            if not i:
                self.layerProperty[self.config.LP['GM']].append(0.0)
                self.layerProperty[self.config.LP['LAPSE']].append(0.0)
            else:
                rho = (amulyr * P) / (properties.R * T)
                dr = abs(zv - self.gas[self.config.C['Z']][i - 1])
                dV = 4.0 * math.pi * (rr**2) * dr
                dM = 1.0e11 * rho * dV
                GdM = self.layerProperty[self.config.LP['GM']][i - 1] + properties.GravConst * dM    # in km3/s2
                self.layerProperty[self.config.LP['GM']].append(GdM)  # mass added as you make way into atmosphere by radius r (times G)
                dT = abs(T - self.gas[self.config.C['T']][i - 1])
                self.layerProperty[self.config.LP['LAPSE']].append(dT / dr)
            # ##set refractivity and index of refraction
            refrlyr = 0.0
            for key in properties.refractivity:
                if key in self.config.C:
                    refrlyr += properties.refractivity[key] * self.gas[self.config.C[key]][i]
            refrlyr = refrlyr * P * (293.0 / T)
            self.layerProperty[self.config.LP['REFR']].append(refrlyr)
            nlyr = refrlyr / 1.0E6 + 1.0
            self.layerProperty[self.config.LP['N']].append(nlyr)

        # ##Now need to normalize GM to planet and calculate scale height (H)
        GMnorm = self.layerProperty[self.config.LP['GM']][iOffset]  # G*(Mass added by p_ref)
        for i, mv in enumerate(self.layerProperty[self.config.LP['GM']]):
            gm = self.config.GM_ref - (mv - GMnorm)
            self.layerProperty[self.config.LP['GM']][i] = gm
            little_g = gm / self.layerProperty[self.config.LP['R']][i]**2
            m_bar = self.layerProperty[self.config.LP['AMU']][i]
            T = self.gas[self.config.C['T']][i]
            self.layerProperty[self.config.LP['H']].append((properties.R * T) / (little_g * m_bar) / 1000.0)
            self.layerProperty[self.config.LP['g']].append(little_g)
        self.layerProperty = np.array(self.layerProperty)

    def computeCloud(self, verbose=False):
        """This computes cloud stuff"""
        print('computeCloud does nothing yet.  This probably wont do anything since gas/cloud/other will get computed in the same tcm')

    def computeGas(self, verbose=False):
        """Computes an atmosphere given stuff"""
        print('computeGas does nothing yet.  This probably wont do anything since gas/cloud/other will get computed in the same tcm')
        print('This will probably just call a external tcm module')

    def __isPresent__(self, c, tiny=1.0E-30):
        """This checks to see if a constituent is there and sets 0.0 or negative values to tiny.  This is generally for log plotting"""
        vsum = 0.0
        present = False
        vnew = []
        for v in c:
            vsum += v
            if v == 0.0:
                vnew.append(tiny)
            elif v < 0.0:
                vnew.append(tiny)
            else:
                vnew.append(v)
                present = True

        return present, vnew
