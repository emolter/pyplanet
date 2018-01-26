# ## This is the file to calculate the radiometric properties of the planets
from __future__ import absolute_import, division, print_function
import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as ss
import sys
import os.path
import prog_path
import utils
import raypath as ray


class Brightness():

    def __init__(self, log=None, verbose=False, plot=False):
        """This calculates the brightness temperature of the planets.
           It must be used with atmosphere and alpha"""
        self.verbose = verbose
        self.plot = plot
        self.log = utils.setupLogFile(log)
        self.layerAlpha = None
        print('\n---Brightness---\n')

    def resetLayers(self):
        self.layerAlpha = None

    def layerAbsorption(self, freqs, atm, alpha):
        self.freqs = freqs
        self.layerAlpha = self.__layerAbsorp__(freqs, atm, alpha)
        if self.plot:
            P = atm.gas[atm.config.C['P']]
            plt.figure('alpha')
            for i, f in enumerate(freqs):
                label = '{:.1f} GHz'.format(f)
                plt.loglog(self.layerAlpha[i], P, label=label)
            v = list(plt.axis())
            v[2] = 100.0 * math.ceil(atm.gas[atm.config.C['P']][-1] / 100.0)
            v[3] = 1.0E-7 * math.ceil(atm.gas[atm.config.C['P']][0] / 1E-7)
            plt.axis(v)
            plt.xlabel(utils.alphaUnit)
            plt.ylabel('P [bars]')
            plt.legend()
            # lgd=plt.legend(loc='upper left',bbox_to_anchor=(1,1))
            # lgd.set_visible(True)  # This is just to remind me...

    def __layerAbsorp__(self, freqs, atm, alpha):
        numLayers = len(atm.gas[0])
        layerAlp = []
        P = atm.gas[atm.config.C['P']]
        T = atm.gas[atm.config.C['T']]
        utils.log(self.log, '{} layers'.format(numLayers), True)
        for layer in range(numLayers):
            if self.verbose:
                print('\r\tAbsorption in layer {}   '.format(layer + 1), end='')
                sys.stdout.flush()
            layerAlp.append(alpha.getAlpha(freqs, T[layer], P[layer], atm.gas[:, layer], atm.config.C, atm.cloud[:, layer],
                            atm.config.Cl, units=utils.alphaUnit))
        layerAlp = np.array(layerAlp).transpose()
        return layerAlp

    def single(self, freqs, atm, b, alpha, orientation=None, taulimit=20.0, discAverage=False, normW4plot=True):
        """This computes the brightness temperature along one ray path"""

        if self.layerAlpha is None:
            self.layerAbsorption(freqs, atm, alpha)
        # get path lengths (ds_layer) vs layer number (num_layer) - currently frequency independent refractivity
        self.path = ray.compute_ds(atm, b, orientation, gtype=None, verbose=self.verbose, plot=self.plot)
        if self.path.ds is None:
            print('Off planet')
            self.Tb = []
            for j in range(len(freqs)):
                self.Tb.append(utils.T_cmb)
            return self.Tb

        # set and initialize arrays
        taus = []
        Tbs = []
        Ws = []
        for j in range(len(freqs)):
            taus.append(0.0)
            Tbs.append(0.0)
            Ws.append(0.0)
        self.tau = [taus]
        self.W = [Ws]
        self.Tb_lyr = [Tbs]

        P = atm.gas[atm.config.C['P']]
        T = atm.gas[atm.config.C['T']]
        z = atm.gas[atm.config.C['Z']]
        self.P = [P[self.path.layer4ds[0]]]
        self.z = [z[self.path.layer4ds[0]]]

        for i in range(len(self.path.ds) - 1):
            ds = self.path.ds[i] * utils.Units[utils.atmLayerUnit] / utils.Units['cm']
            taus = []
            Ws = []
            Tbs = []
            ii = self.path.layer4ds[i]
            ii1 = self.path.layer4ds[i + 1]
            T1 = T[ii1]
            T0 = T[ii]
            self.P.append((P[ii] + P[ii1]) / 2.0)
            self.z.append((z[ii] + z[ii1]) / 2.0)

            if self.layerAlpha is None:
                print("is None at ", i)
            for j, f in enumerate(freqs):
                if not alpha.config.Doppler:
                    a1 = self.layerAlpha[j][ii1]
                    a0 = self.layerAlpha[j][ii]
                else:
                    fshifted = [[f / self.path.doppler[i]], [f / self.path.doppler[i + 1]]]
                    print('\rdoppler corrected frequency at layer', i, end='')
                    a1 = alpha.getAlpha(fshifted[0], T[ii1], P[ii1], atm.gas[:, ii1], atm.config.C, atm.cloud[:, ii1],
                                        atm.config.Cl, units=utils.alphaUnit)
                    a0 = alpha.getAlpha(fshifted[1], T[ii], P[ii], atm.gas[:, ii], atm.config.C, atm.cloud[:, ii],
                                        atm.config.Cl, units=utils.alphaUnit)
                dtau = (a0 + a1) * ds / 2.0
                taus.append(self.tau[i][j] + dtau)         # this is tau_(i+1)
                if discAverage is True:
                    Ws.append(2.0 * a1 * ss.expn(2, taus[j]))  # this is W_(i+1) for disc average
                else:
                    Ws.append(a1 * math.exp(-taus[j]))  # this is W_(i+1) for non disc average
                dTb = (T1 * Ws[j] / scriptR(T1, freqs[j]) + T0 * self.W[i][j] / scriptR(T0, freqs[j])) * ds / 2.0
                Tbs.append(self.Tb_lyr[i][j] + dTb)
            self.tau.append(taus)
            self.W.append(Ws)
            self.Tb_lyr.append(Tbs)

        # final spectrum
        self.Tb = []
        for j in range(len(freqs)):
            top_Tb_lyr = self.Tb_lyr[-1][j]
            if top_Tb_lyr < utils.T_cmb:
                top_Tb_lyr = utils.T_cmb
            self.Tb.append(top_Tb_lyr)
        self.tau = np.array(self.tau).transpose()
        self.W = np.array(self.W).transpose()
        self.Tb_lyr = np.array(self.Tb_lyr).transpose()
        self.P = np.array(self.P)
        self.z = np.array(self.z)

        if self.plot:
            # ####-----Weigthing functions
            plt.figure('radtran')
            plt.subplot(121)
            for i, f in enumerate(freqs):
                # label=r'$\tau$: %.1f GHz' % (f)
                # plt.semilogy(self.tau[i],self.P,label=label)
                if normW4plot:
                    wplot = self.W[i] / np.max(self.W[i])
                else:
                    wplot = self.W[i]
                label = (r'$W$: {:.1f} GHz').format(f)
                label = (r'{:.1f} cm').format(30.0 / f)
                if len(wplot) == len(self.P):
                    plt.semilogy(wplot, self.P, label=label, linewidth=3)
                else:
                    print("Not plotted since wplot is length {} and P is length {}".format(len(wplot), len(self.P)))
            plt.legend()
            plt.axis(ymin=100.0 * math.ceil(np.max(self.P) / 100.0), ymax=1.0E-7 * math.ceil(np.min(self.P) / 1E-7))
            plt.ylabel('P [bars]')
            # ####-----Alpha
            plt.figure('alpha')
            for i, f in enumerate(freqs):
                label = (r'$\alpha$: {:.1f} GHz').format(f)
                label = (r'{:.1f} cm').format(30.0 / f)
                pl = list(self.layerAlpha[i])
                del pl[0]
                # delete because alpha is at the layer boundaries, so there are n+1 of them
                if len(pl) == len(self.P):
                    plt.loglog(pl, self.P, label=label)
                else:
                    print("Not plotted since wplot is length {} and P is length {}".format(len(pl), len(self.P)))
            plt.legend()
            v = list(plt.axis())
            v[2] = 100.0 * math.ceil(np.max(self.P) / 100.0)
            v[3] = 1.0E-7 * math.ceil(np.min(self.P) / 1E-7)
            plt.axis(v)
            plt.ylabel('P [bars]')
            # ####-----Brightness temperature
            plt.figure('brightness')
            lt = '-'
            if (len(self.Tb) == 1):
                lt = 'o'
            plt.plot(freqs, self.Tb, lt)
            plt.xlabel('Frequency [GHz]')
            plt.ylabel('Brightness temperature [K]')

        del taus, Tbs, Ws
        return self.Tb

    def savertm(self, tag=None, path='Output'):
        if tag is None:
            filename = None
        else:
            filename = 'alpha_' + tag + '.out'
        self.saveAlpha(filename, path)
        if tag is None:
            filename = None
        else:
            filename = 'wgt_' + tag + '.out'
        self.saveWeight(filename, path)
        if tag is None:
            filename = None
        else:
            filename = 'tau_' + tag + '.out'
        self.saveTau(filename, path)
        if tag is None:
            filename = None
        else:
            filename = 'tblayer_' + tag + '.out'
        self.saveTblayer(filename, path)

    def saveAlpha(self, filename=None, path='.'):
        if filename is None:
            filename = 'alpha.out'
        os.path.join(path, filename)
        fp = open(filename, 'w')
        s = '#P  \tz  \t'
        for f in self.freqs:
            s += '{:.2f}\t'.format(f)
        s += 'GHz\n'
        fp.write(s)
        for j in range(len(self.P)):
            s = ('{}\t{:.2f}\t').format(repr(self.P[j]), self.z[j])
            for i in range(len(self.freqs)):
                s += '{}\t'.format(repr(self.layerAlpha[i][j]))
            s += '\n'
            fp.write(s)
        s = ('%s (%d x %d)').format(filename, i + 1, j + 1)

    def saveWeight(self, norm=False, filename=None, path='.'):
        if filename is None:
            filename = 'wgt.out'
        fp = open(filename, 'w')
        s = '#P  \tz  \t'
        for f in self.freqs:
            s += ('{:.2f}\t').format(f)
        s = s.strip() + 'GHz\n'
        fp.write(s)
        scale = []
        for i in range(len(self.freqs)):
            if norm:
                scale.append(np.max(self.W[i]))
            else:
                scale.append(1.0)
        for j in range(len(self.P)):
            s = ('%s\t%.2f\t').format(repr(self.P[j]), self.z[j])
            for i in range(len(self.freqs)):
                s += ('%s\t').format(repr(self.W[i][j] / scale[i]))
            s = s.strip() + '\n'
            fp.write(s)
        s = ('%s (%d x %d)').format(filename, i + 1, j + 1)
        return s

    def saveTau(self, filename=None, path='.'):
        if filename is None:
            filename = 'tau.out'
        os.path.join(path, filename)
        fp = open(filename, 'w')
        s = '#P  \tz  \t'
        for f in self.freqs:
            s += '{:.2f}\t'.format(f)
        s += 'GHz\n'
        fp.write(s)
        for j in range(len(self.P)):
            s = ('%s\t%.2f\t').format(repr(self.P[j]), self.z[j])
            for i in range(len(self.freqs)):
                s += ('%s\t').format(repr(self.tau[j]))
            s += '\n'
            fp.write(s)
        s = ('%s (%d x %d)').format(filename, i + 1, j + 1)
        return s

    def saveTblayer(self, filename=None, path='.'):
        if filename is None:
            filename = 'tblayer.out'
        os.path.join(path, filename)
        fp = open(filename, 'w')
        s = '#P  \tz  \t'
        for f in self.freqs:
            s += ('%.2f\t').format(f)
        s += 'GHz\n'
        fp.write(s)
        for j in range(len(self.P)):
            s = ('{}\t{:.2f}\t').format(repr(self.P[j]), self.z[j])
            for i in range(len(self.freqs)):
                s += ('%s\t').format(repr(self.Tb_lyr[j]))
            s += '\n'
            fp.write(s)
        s = ('%s (%d x %d)').format(filename, i + 1, j + 1)
        return s


def scriptR(T, freq):
    """See Janssen pg 7"""
    a = (utils.hP * freq * utils.Units[utils.processingFreqUnit]) / (utils.kB * T)
    R = (math.exp(a) - 1.0) / a
    return R
