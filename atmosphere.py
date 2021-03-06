from __future__ import absolute_import, division, print_function
import math
import matplotlib.pyplot as plt
import numpy as np
import utils
import sys
import os

# ##local imports
import prog_path
import config as pcfg
import raypath as ray
import chemistry
import regrid
import state_variables


class Atmosphere:
    def __init__(self, planet, mode='normal', config='config.par', log=None, **kwargs):
        """reads/computes atmospheres.  This returns:
               self.gas
               self.cloud
               self.layerProperty
            on the appropriate grid."""

        self.planet = planet.capitalize()
        kwargs = state_variables.init_state_variables(mode, **kwargs)
        self.state_vars = kwargs.keys()
        self.set_state(set_mode='init', **kwargs)
        if self.verbose:
            print('\n---Atmosphere of {}---'.format(planet))
        self.logFile = utils.setupLogFile(log)

        if type(config) == str:
            config = os.path.join(self.planet, config)
            config = pcfg.planetConfig(self.planet, configFile=config, log=log)
        self.config = config

        # ##Create function dictionaries
        self.gasGen = {}
        self.gasGen['read'] = self.readGas
        self.cloudGen = {}
        self.cloudGen['read'] = self.readCloud
        self.propGen = {}
        self.propGen['compute'] = self.computeProp

        if self.verbose == 'loud':
            print('Planet ' + self.planet)
            self.config.display()
        if self.config.gasType == 'read':  # this assumes that cloudType is then also 'read'
            utils.log(self.logFile, '\tReading from: ' + self.config.filename, self.verbose)
            utils.log(self.logFile, '\tAtmosphere file:  ' + self.config.gasFile, self.verbose)
            utils.log(self.logFile, '\tCloud file:  ' + self.config.cloudFile, self.verbose)

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

        # ## Generate gas profile
        if gasType not in self.gasGen.keys():
            print('Error:  No such gasType: ', gasType)
            return 0
        else:
            self.gasGen[gasType]()
        self.nAtm = len(self.gas[0])

        if self.batch_mode:
            return self.nAtm

        self.chem = {}
        for key in self.config.C:
            if key in ['P', 'T', 'Z', 'DZ']:
                continue
            self.chem[key] = chemistry.ConstituentProperties(key)

        # ## Generate cloud profile
        if cloudType not in self.cloudGen.keys():
            print('Error:  No such cloudType: ', cloudType)
            return 0
        else:
            self.cloudGen[cloudType]()

        # ## Put onto common grid
        if self.verbose:
            print("Regrid:  {}".format(regridType))
        regridded = regrid.regrid(self, regridType=regridType, Pmin=Pmin, Pmax=Pmax)
        self.nAtm = len(self.gas[0])

        if tweak:  # This loads and calls the module as given in the config.par tweakmodule parameter
            self.tweakAtm()

        # ## Compute other parameters that are needed
        if self.config.otherType not in self.propGen.keys():
            print('Error:  no such otherTpe: ', otherType)
            return 0
        else:
            self.propGen[otherType]()

        angularDiameter = 2.0 * math.atan(self.layerProperty[self.config.LP['R']][0] / self.config.distance)
        if self.verbose == 'loud':
            print('angular radius = {} arcsec'.format(utils.r2asec(angularDiameter / 2.0)))

        # ## Plot data
        if self.plot:
            self.plotTP()
            self.plotGas()
            self.plotCloud()
            self.plotProp()

        return self.nAtm

    def set_state(self, set_mode='set', **kwargs):
        """
        set_mode:  'set' or 'init', if set, checks list
        """
        for k, v in kwargs.iteritems():
            if isinstance(v, str):
                v = v.lower()
            if k in self.state_vars:
                setattr(self, k, v)
                if set_mode == 'set':
                    print('Setting {} to {}'.format(k, v))
            else:
                if set_mode == 'set':
                    print('state_var [{}] not found.'.format(k))
        if set_mode == 'init' and self.verbose == 'loud':
            self.show_state()

    def show_state(self):
        print("Atmosphere state variables")
        for k in self.state_vars:
            print('\t{}:  {}'.format(k, getattr(self, k)))

    def readGas(self, gasFile=None, Cdict=None):
        """Reads gas profile file as self.gas"""

        if gasFile is None:
            gasFile = self.config.gasFile
        if Cdict is None:
            Cdict = self.config.C
        if gasFile == 'batch':
            self.batch_mode = True
            return 0
        else:
            self.batch_mode = False
        gasFile = os.path.join(self.config.path, gasFile)

        if self.verbose:
            print('Reading constituents from {}'.format(gasFile))
            print('\tUsing atmsopheric component:  ', end='')
        self.gas = []
        for k in Cdict:
            self.gas.append([])
            if self.verbose:
                print(k + '  ', end='')
        self.nConstituent = len(Cdict)
        try:
            fp = open(gasFile, "r")
        except IOError:
            print(gasFile + ' was not found - returning no gas profile\n\n')
            raise IOError
        if self.verbose:
            print(' ')
        expected_number = utils.get_expected_number_of_entries(fp)
        for line in fp:
            cval = utils.get_data_from(line)
            if cval is None or len(cval) != expected_number:
                continue
            for n in range(self.nConstituent):
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

        if self.verbose:
            print('\tRead ' + str(self.nGas) + ' lines')
        fp.close()
        self.gas = np.array(self.gas)
        # ## Check that P is monotonically increasing
        monotonic = np.all(np.diff(self.gas[self.config.C['P']]) > 0.0)
        if not monotonic:
            self.gas = np.fliplr(self.gas)
            monotonic = np.all(np.diff(self.gas[self.config.C['P']]) > 0.0)
        if not monotonic:
            print("Error in ", gasFile)
            raise ValueError("Pressure not monotonically increasing.")

        # ## Renormalize so that deepest z is 0 and set DZ
        self.renorm_z('gas')

        return self.nGas

    def readCloud(self, cloudFile=None, Cldict=None):
        """Reads in cloud data if we have it..."""

        if cloudFile is None:
            cloudFile = self.config.cloudFile
        if Cldict is None:
            Cldict = self.config.Cl
        if cloudFile == 'batch':
            self.batch_mode = True
            return 0
        else:
            self.batch_mode = False
        cloudFile = os.path.join(self.config.path, cloudFile)

        if self.verbose:
            print('Reading clouds from {}'.format(cloudFile))
            print('\tUsing cloud components:  ', end='')
        self.cloud = []
        for k in Cldict:
            self.cloud.append([])
            if self.verbose:
                print(k, '  ', end='')
        self.nParticulate = len(Cldict.keys())
        try:
            fp = open(cloudFile, "r")
        except IOError:
            print(cloudFile, ' was not found - returning no clouds\n\n')
            raise IOError
        if self.verbose:
            print(' ')
        expected_number = utils.get_expected_number_of_entries(fp)
        for line in fp:
            cval = utils.get_data_from(line)
            if cval is None or len(cval) != expected_number:
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

        if self.verbose:
            print('\tRead ', str(self.nCloud), ' lines')
        fp.close()
        self.cloud = np.array(self.cloud)
        # ## Check that P is monotonically increasing
        monotonic = np.all(np.diff(self.cloud[self.config.Cl['P']]) > 0.0)
        if not monotonic:
            self.cloud = np.fliplr(self.cloud)
            monotonic = np.all(np.diff(self.cloud[self.config.Cl['P']]) > 0.0)
        if not monotonic:
            print("Error in ", cloudFile)
            raise ValueError("Pressure not monotonically increasing")

        # ## Renormalize so that deepest z is 0 and set DZ
        self.renorm_z('cloud')

        return self.nCloud

    def renorm_z(self, gctype):
        if gctype == 'gas':
            ind = self.config.C
            atm_gc = self.gas
        else:
            ind = self.config.Cl
            atm_gc = self.cloud
        zDeep = atm_gc[ind['Z']][-1]
        for i in range(len(atm_gc[ind['Z']])):
            atm_gc[ind['Z']][i] -= zDeep

        # put in DZ
        dz = np.abs(np.diff(atm_gc[ind['Z']])) * 1.0E5  # convert from km to cm (so no unit checking!!!)
        atm_gc[ind['DZ']] = np.append(np.array([0.0]), dz)

    def tweakAtm(self):
        """Tweaks the atmosphere data..."""
        if self.batch_mode:
            return 0

        # Import tweakmodule
        sys.path.append(self.config.path)
        try:
            __import__(self.config.tweakmodule)
            tweakModule = sys.modules[self.config.tweakmodule]
        except SyntaxError:
            utils.log(self.logFile, "Syntax Error:  check " + self.config.tweakmodule, True)
            raise ValueError("Error in tweakmodule")

        # Run module then log
        self.tweakComment, self.gas, self.cloud = tweakModule.modify(self.gas, self.cloud, self.config.C, self.config.Cl)
        if self.verbose:
            print('---tweakComment')
            print(self.tweakComment)
            print('---')
        utils.log(self.logFile, self.tweakComment, False)
        tf = os.path.join(self.config.path, self.config.tweakmodule + '.py')
        tp = open(tf, 'r')
        dt = tp.read()
        utils.log(self.logFile, '======================' + tf + '=====================', False)
        utils.log(self.logFile, dt, False)
        utils.log(self.logFile, '====================================================================', False)
        tp.close()

    def scaleAtm(self, scale_info='Scratch/scale.dat', plot_diff=True):
        """
        This is a built-in tweak module.
        """
        if isinstance(scale_info, str):
            import alpha
            col, scale_info = alpha.read_scalefile(scale_info)
        else:
            col = scale_info.keys()

        if len(scale_info[col[0]]) != self.nAtm:
            print("Warning - scale file doesn't match atmosphere.  Not applying.")
            return None

        if plot_diff:
            color_seq = ['b', 'k', 'r', 'm', 'c']
            clr = {}
            plt.figure('Scale difference')
            for i, gas in enumerate(col):
                if gas.lower() != 'p':
                    clr[gas] = color_seq[i]
                    present, g = self.is_present(self.gas[self.config.C[gas.upper()]])
                    plt.loglog(g, self.gas[self.config.C['P']], color=clr[gas], linestyle='--', label=gas + ' before')

        for i in range(self.nAtm):
            for gas in col:
                if gas.lower() != 'p':
                    self.gas[self.config.C[gas.upper()]][i] *= scale_info[gas][i]

        if plot_diff:
            for i, gas in enumerate(col):
                if gas.lower() != 'p':
                    present, g = self.is_present(self.gas[self.config.C[gas.upper()]])
                    plt.loglog(g, self.gas[self.config.C['P']], color=clr[gas], linestyle='-', label=gas + ' after')
            self.frame_plot('Fractional Abundance')

    def computeProp(self):
        """This module computes derived atmospheric properties (makes self.layerProperty)"""
        if self.batch_mode:
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
        if self.verbose == 'loud':
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
            for key in self.chem:
                amulyr += self.chem[key].amu * self.gas[self.config.C[key]][i]
            self.layerProperty[self.config.LP['AMU']].append(amulyr)
            # ##set GM pre-calc (normalized further down) and get lapse rate
            if not i:
                self.layerProperty[self.config.LP['GM']].append(0.0)
                self.layerProperty[self.config.LP['LAPSE']].append(0.0)
                self.layerProperty[self.config.LP['LAPSEP']].append(0.0)
            else:
                rho = (amulyr * P) / (chemistry.R * T)
                dr = abs(zv - self.gas[self.config.C['Z']][i - 1])
                dV = 4.0 * math.pi * (rr**2) * dr
                dM = 1.0e11 * rho * dV
                GdM = self.layerProperty[self.config.LP['GM']][i - 1] + chemistry.GravConst * dM    # in km3/s2
                self.layerProperty[self.config.LP['GM']].append(GdM)  # mass added as you make way into atmosphere by radius r (times G)
                dT = abs(T - self.gas[self.config.C['T']][i - 1])
                dP = abs(P - self.gas[self.config.C['P']][i - 1])
                self.layerProperty[self.config.LP['LAPSE']].append(dT / dr)
                self.layerProperty[self.config.LP['LAPSEP']].append(dT / dP)
            # ##set refractivity and index of refraction
            refrlyr = 0.0
            for key in self.chem:
                refrlyr += self.chem[key].refractivity(T=T) * self.gas[self.config.C[key]][i]
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
            self.layerProperty[self.config.LP['H']].append((chemistry.R * T) / (little_g * m_bar) / 1000.0)
            self.layerProperty[self.config.LP['g']].append(little_g)
        self.layerProperty = np.array(self.layerProperty)

    def is_present(self, c, tiny=1.0E-30):
        """This checks to see if a constituent is there and sets 0.0 or negative values to tiny.
           This is generally for log plotting."""
        v = [tiny if x <= tiny else x for x in c]
        present = bool(len(np.where(np.array(v) > tiny)[0]))
        return present, v

    def frame_plot(self, xlabel, show_legend=True):
        v = list(plt.axis())
        if v[0] < 1E-10:
            v[0] = 1E-10
        v[2] = 100.0 * math.ceil(self.gas[self.config.C['P']][-1] / 100.0)
        v[3] = 1.0E-7 * math.ceil(self.gas[self.config.C['P']][0] / 1E-7)
        plt.axis(v)
        plt.ylabel('P [bars]')
        plt.xlabel(xlabel)
        if show_legend:
            plt.legend()

    def plotTP(self, plot='auto'):
        """Plot the T-P profile"""
        if plot == 'auto':
            plt.figure(self.planet + ': T-P')
        plt.title(self.planet + ':  T-P profile')
        plt.loglog(self.gas[self.config.C['T']], self.gas[self.config.C['P']], label='T')
        self.frame_plot('T [K]', show_legend=False)

    def plotCloud(self, dontPlot=['Z', 'P', 'T', 'DZ'], plot='auto'):
        """Plots the clouds"""
        if plot == 'auto':
            plt.figure(self.planet + ': clouds')
        plt.title(self.planet + ': clouds')
        for cloud in self.config.Cl:
            present, cl = self.is_present(self.cloud[self.config.Cl[cloud]])
            if cloud in dontPlot or not present:
                continue
            plt.loglog(cl, self.cloud[self.config.Cl['P']], label=cloud)
        self.frame_plot(r'Density [g/cm$^3$]')

    def plotGas(self, dontPlot=['Z', 'P', 'T', 'DZ'], plot='auto'):
        """Plots the constituents"""
        if plot == 'auto':
            plt.figure(self.planet + ': gas')
        plt.title(self.planet + ': gas')
        for gas in self.config.C:
            present, g = self.is_present(self.gas[self.config.C[gas]])
            if gas in dontPlot or not present:
                continue
            plt.loglog(g, self.gas[self.config.C['P']], label=gas)
        self.frame_plot('Fractional Abundance')

    def plotProp(self, dontPlot=['Z', 'P', 'T'], plot='auto'):
        if plot == 'auto':
            plt.figure(self.planet + ': Properties')
        plt.title(self.planet + ': other')
        for other in self.config.LP:
            present, g = self.is_present(self.layerProperty[self.config.LP[other]])
            if other in dontPlot or not present:
                continue
            plt.loglog(g, self.gas[self.config.C['P']], label=other)
        self.frame_plot('Property Value')
