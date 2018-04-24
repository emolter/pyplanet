#  ## This is the class to calculate the microwave properties of the constituents
from __future__ import absolute_import, division, print_function
import os
import sys
import numpy as np
import prog_path
import utils
import config as pcfg
import state_variables


class Alpha:
    def __init__(self, config=None, log=None, **kwargs):
        """Reads in absorption formalisms
           Note that they are all in GHz"""

        kwargs = state_variables.init_state_variables('normal', **kwargs)
        self.state_vars = kwargs.keys()
        self.set_state(set_mode='init', **kwargs)
        if self.verbose:
            self.show_state('Alpha')
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
        if self.use_existing_alpha:
            self.read_existing_alpha()
        else:
            self.formalisms()
        if self.generate_alpha:
            np.savez('Scratch/constituents', alpha_dict=self.config.constituent_alpha, alpha_sort=self.ordered_constituents)
            self.fp_gen_alpha = open('Scratch/absorb.dat', 'w')

    def complete_generate_alpha(self, n_freq, n_layer):
        self.fp_gen_alpha.close()
        data = []
        with open('Scratch/absorb.dat', 'r') as fp:
            for i, line in enumerate(fp):
                if not i % n_freq:
                    if i:
                        data.append(layer)
                    layer = []
                v = [float(x) for x in line.split()]
                layer.append(v)
            data.append(layer)
        data = np.array(data)
        np.save('Scratch/absorb', data)
        os.remove('Scratch/absorb.dat')

    def read_existing_alpha(self):
        self.alpha_data = np.load('Scratch/absorb.npy')

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
        self.ordered_constituents = sorted(self.constituent.keys())
        utils.log(self.log, 'Using modules:', True)
        for k in self.constituent:
            utils.log(self.log, '\t' + k + ':  ' + self.constituent[k], True)

    def getAlpha(self, freqs, layer, atm, units='invcm', plot=None):
        """This is a wrapper to get the absorption coefficient, either from calculating from formalisms
           or reading from file"""
        if self.use_existing_alpha:
            return self.get_alpha_from_file(freqs, layer, units, plot)
        else:
            P = atm.gas[atm.config.C['P']][layer]
            T = atm.gas[atm.config.C['T']][layer]
            gas = atm.gas[:, layer]
            cloud = atm.cloud[:, layer]
            return self.get_alpha_from_calc(freqs, T, P, gas, atm.config.C, cloud, atm.config.Cl, units, plot)

    def get_alpha_from_file(self, freqs, layer, units='invcm', plot=None):
        totalAbsorption = np.zeros_like(freqs)
        for i in range(len(freqs)):
            totalAbsorption[i] = self.alpha_data[layer, i, -1]
        return totalAbsorption

    def get_alpha_from_calc(self, freqs, T, P, gas, gas_dict, cloud, cloud_dict, units='invcm', plot=None):
        """This gets the total absoprtion coefficient from gas.  It assumes the correct frequency units, but maybe should correct that.
           Returns total absorption at that layer."""
        absorb = []
        for k in self.ordered_constituents:
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
        if self.generate_alpha:
            self.write_layer(absorb, totalAbsorption)
        return totalAbsorption

    def write_layer(self, absorb, totalAbsorption):
        for i in range(len(absorb)):
            s = ''
            for a2 in absorb[i]:
                s += '{} '.format(a2)
            s += '{}\n'.format(totalAbsorption[i])
            self.fp_gen_alpha.write(s)

    def set_state(self, set_mode='set', **kwargs):
        for k, v in kwargs.iteritems():
            if k in self.state_vars:
                setattr(self, k, v)
                if set_mode == 'set':
                    print('Setting {} to {}'.format(k, v))
            else:
                if set_mode == 'set':
                    print('state_var [{}] not found.'.format(k))

    def show_state(self, stype):
        print("{} state variables".format(stype))
        for k in self.state_vars:
            print('\t{}:  {}'.format(k, getattr(self, k)))
