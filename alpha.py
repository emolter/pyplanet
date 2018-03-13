#  ## This is the class to calculate the microwave properties of the constituents
from __future__ import absolute_import, division, print_function
import math
import string
import matplotlib.pyplot as plt
import os
import os.path
import sys
import numpy as np
import prog_path
import utils
import config as pcfg


class Alpha:
    def __init__(self, config=None, log=None, verbose=False, plot=False):
        """Reads in absorption formalisms
           Note that they are all in GHz"""

        self.verbose = verbose
        self.plot = plot
        self.log = utils.setupLogFile(log)

        # Get possible constituents
        possible = []
        self.constituentsAreAt = prog_path.constituentPath
        s = 'Reading in absorption modules from ' + self.constituentsAreAt + '\n'
        utils.log(self.log, s, verbose)
        for d in os.listdir(self.constituentsAreAt):
            fnd = os.path.join(self.constituentsAreAt, d)
            if os.path.isdir(fnd):
                possible.append(d)
        # Import used ones - note this dynamically imports the absorption modules.  It checks that directory's use.txt file.
        self.constituent = {}
        self.absorptionModule = {}
        for c in possible:
            fn = os.path.join(self.constituentsAreAt, c, 'use.txt')
            try:
                fp = open(fn, 'r')
            except IOError:
                utils.log(self.log, 'No file ' + fn, True)
                continue
            absorber = fp.readline().strip()
            testabs = absorber.split('.')
            if len(testabs) == 2:
                absorber = testabs[0]
            fp.close()
            constituentPath = os.path.join(self.constituentsAreAt, c)
            if string.lower(absorber) != 'none':
                sys.path.append(constituentPath)
                try:
                    __import__(absorber)
                    self.absorptionModule[c] = sys.modules[absorber]
                    self.constituent[c] = absorber
                except ImportError:
                    s = "WARNING:  CAN'T LOAD " + absorber + '\n'
                    print(s * 3)
                    utils.log(self.log, "Can't load " + absorber, True)
        utils.log(self.log, 'Using modules:', True)
        for k in self.constituent:
            utils.log(self.log, '\t' + k + ':  ' + self.constituent[k], True)

        # get config
        if type(config) == str:
            config = pcfg.planetConfig(self.planet, configFile=config, log=log)
        self.config = config

        # copy config back into otherPar
        self.otherPar = {}
        self.otherPar['h2state'] = self.config.h2state
        self.otherPar['h2newset'] = self.config.h2newset
        self.otherPar['water'] = self.config.water_p
        self.otherPar['ice'] = self.config.ice_p
        self.otherPar['nh4sh'] = self.config.nh4sh_p
        self.otherPar['nh3ice'] = self.config.nh3ice_p
        self.otherPar['h2sice'] = self.config.h2sice_p
        self.otherPar['ch4'] = self.config.ch4_p

    def test(self, f=[1., 2., 3., 4., 5.], T=300.0, P=1.0, X=[0.9, 0.1, 0.001], D=None,
             otherPar=None, units='dBperkm', fignum=10, plot=True):
        """This just tests things - as well as reminds me how it works"""
        if plot is None:
            plot = self.plot
        if otherPar is None:
            otherPar = {'h2state': 'e', 'h2_newset': True}
        print('Testing alpha modules:')
        print('\tT = ', T)
        print('\tP = ', P)
        if f == 0:
            freq = []
            for i in range(500):
                freq.append(float(i))
            f = freq
        if D is None:
            D = {}
            for k in self.constituent:
                D[string.upper(k)] = 2
            D['H2'] = 0
            D['HE'] = 1
            D['CH4'] = 2
        i = 0
        self.a = []
        if plot:
            plt.figure(fignum)
        for k in self.constituent:
            print('----------------------' + k + ':  ' + self.constituent[k])
            path = os.path.join(self.constituentsAreAt, k)
            self.a.append(self.absorptionModule[k].alpha(f, T, P, X, D, otherPar, units=units, path=path, verbose=True))
            if plot:
                plt.semilogy(f, self.a[i], label=k)
            i += 1
        if plot:
            plt.xlabel('Frequency [GHz]')
            plt.ylabel('Absorption [dB/km]')
            plt.legend()
        print('Done test.  See plot {:d} -------'.format(fignum))

    def getAlpha(self, freqs, T, P, gas, gas_dict, cloud, cloud_dict, units='invcm', plot=None):
        """This gets the total absoprtion coefficient from gas.  It assumes the correct frequency units, but maybe should correct that."""
        absorb = []
        for k in self.constituent:
            path = os.path.join(self.constituentsAreAt, k)
            if k[0:4].lower() == 'clou':
                X = cloud
                D = cloud_dict
            else:
                X = gas
                D = gas_dict
            absorb.append(self.absorptionModule[k].alpha(freqs, T, P, X, D, self.otherPar, units=units, path=path, verbose=self.verbose))
        absorb = np.array(absorb)
        absorb = absorb.transpose()
        totalAbsorption = np.zeros_like(freqs)
        for i in range(len(freqs)):
            totalAbsorption[i] = absorb[i].sum()
        return totalAbsorption
