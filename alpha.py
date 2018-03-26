#  ## This is the class to calculate the microwave properties of the constituents
from __future__ import absolute_import, division, print_function
import os
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
        self.formalisms()

    def formalisms(self):
        # Get possible constituents
        self.constituentsAreAt = prog_path.constituentPath
        s = 'Reading in absorption modules from ' + self.constituentsAreAt + '\n'
        utils.log(self.log, s, self.verbose)
        # Import used ones - note this dynamically imports the absorption modules.
        self.constituent = {}
        self.absorptionModule = {}
        for c in self.config.constituent_alpha:
            absorber = self.config.constituent_alpha[c]
            if absorber is None:
                continue
            constituentPath = os.path.join(self.constituentsAreAt, c)
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
