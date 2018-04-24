#  This is the 'executive' class for planets
from __future__ import absolute_import, division, print_function
import math
import matplotlib.pyplot as plt
import numpy as np
import datetime
import prog_path
import atmosphere as atm
import config as pcfg
import alpha
import brightness as bright
import data_handling
import utils
import fileIO
import state_variables
import os

version = '2.3'


class Planet:
    def __init__(self, name, mode='normal', config='planet', **kwargs):
        """This is the 'executive function class to compute overall planetary emission.
           For both mode and kwargs look at state_variables.py
           Inputs:
                name:  'Jupiter', 'Saturn', 'Uranus', 'Neptune'
                config:  config file name.  If 'planet' sets to <name>/config.par
                mode:  '[normal]/batch/mcmc/generate_alpha/use_alpha'  sets up for various special modes
                kwargs: 'verbose' and 'plot' (and other state_vars - see show_state())"""

        planetList = ['Jupiter', 'Saturn', 'Neptune', 'Uranus']
        self.planet = name.capitalize()
        self.header = {}
        self.imrow = False
        self.freqs = None
        self.freqUnit = None
        self.b = None
        self.imSize = None

        print('Planetary modeling  (ver {})'.format(version))

        if self.planet not in planetList:
            print("{} not found.".format(self.planet))
            return

        # Set up state_variables
        kwargs = state_variables.init_state_variables(mode.lower(), **kwargs)
        self.state_vars = kwargs.keys()
        self.set_state(set_mode='init', **kwargs)
        if not self.super_quiet:
            self.show_state('Planet')

        #  ##Set up log file
        if self.write_log_file:
            runStart = datetime.datetime.now()
            self.logFile = 'Logs/{}_{}.log'.format(self.planet, runStart.strftime("%Y%m%d_%H%M"))
            self.log = utils.setupLogFile(self.logFile, not self.super_quiet)
            utils.log(self.log, self.planet + ' start ' + str(runStart), not self.super_quiet)
        else:
            self.log = None

        #  ## Get config
        if config.lower() == 'planet':
            config = self.planet + '/config.par'
        if not self.super_quiet:
            print('Reading config file:  ', config)
            print("\t'print(x.config.show())' to see config parameters (where x in the instance name, e.g. 'j').")
        self.config = pcfg.planetConfig(self.planet, configFile=config, log=self.log)

        #  ## Create atmosphere:  attributes are self.atm.gas, self.atm.cloud and self.atm.layerProperty
        self.atm = atm.Atmosphere(self.planet, config=self.config, log=self.log, **kwargs)
        self.atm.run()

        #  ## Read in absorption modules:  to change absorption, edit files under /constituents'
        self.alpha = alpha.Alpha(config=self.config, log=self.log, **kwargs)

        #  ## Next compute radiometric properties - initialize bright and return data class
        self.bright = bright.Brightness(log=self.log, **kwargs)
        self.data_return = data_handling.DataReturn()

        # ## Create fileIO class
        self.fIO = fileIO.FileIO(self.output_type)

    def run(self, freqs='reuse', b=[0.0, 0.0], freqUnit='GHz', block=[1, 1], orientation=None):
        """Runs the model to produce the brightness temperature, weighting functions etc etc
            freqs:  frequency request as set in set_freq.  If 'reuse' it won't recompute absorption/layer (allows many b)
            b:  "impact parameter" request as set in set_b
            freqUnit:  unit that freqs is in
            block:  blocks to produce image (related to memory error...)
            orientation:  orientation vector of planet"""

        #  ##Set freqs
        if self.use_existing_alpha or self.scale_existing_alpha:
            freqs_read = np.load('Scratch/freqs.npy')
            freqs = [f for f in freqs_read]
            if not self.super_quiet:
                print("Setting frequencies to ", freqs)
        reuse = False
        if isinstance(freqs, str) and freqs == 'reuse':
            if self.freqs is None:
                raise ValueError('Must set frequencies.')
            reuse = True
            freqs = self.freqs
            freqUnit = self.freqUnit
        else:
            freqs, freqUnit = self.set_freq(freqs, freqUnit)
            self.bright.resetLayers()
        self.data_return.f = freqs

        #  ##Set b, etc
        b, outType = self.set_b(b, block)
        self.data_return.b = b
        if outType == 'Image' and len(freqs) > 1:
            print('Warning:  Image must be at only one frequency')
            print('Using {} {}'.format(freqs[0], freqUnit))
            self.freqs = list(freqs[0])
            freqs = self.freqs
        if self.verbose:
            print('outType = {}'.format(outType))

        #  ##Start
        runStart = datetime.datetime.now()
        self.Tb = []
        btmp = ''
        self.rNorm = None
        self.tip = None
        self.rotate = None
        if outType == 'Image':  # We now treat it as an image at one frequency
            if self.verbose:
                print('imgSize = {} x {}'.format(self.imSize[0], self.imSize[1]))
            imtmp = []
            if abs(block[1]) > 1:
                btmp = '_{:02d}of{:02d}'.format(block[0], abs(block[1]))
            else:
                btmp = ''

        for i, bv in enumerate(b):
            if self.verbose:
                print('{} of {} (view [{:.4f}, {:.4f}])  '.format(i + 1, len(b), bv[0], bv[1]), end='')
            Tbt = self.bright.single(freqs, self.atm, bv, self.alpha, orientation, discAverage=(self.bType == 'disc'))
            if self.bright.travel is not None:
                if self.rNorm is None:
                    self.rNorm = self.bright.travel.rNorm
                if self.tip is None:
                    self.tip = self.bright.travel.tip
                if self.rotate is None:
                    self.rotate = self.bright.travel.rotate
            if outType == 'Image':
                imtmp.append(Tbt[0])
                if not (i + 1) % self.imSize[0]:
                    self.Tb.append(imtmp)
                    imtmp = []
            else:
                self.Tb.append(Tbt)
        self.data_return.Tb = self.Tb
        self.data_return.header = self.header
        missed_planet = self.rNorm is None

        if self.generate_alpha:
            self.alpha.complete_generate_alpha(len(freqs), self.atm.nAtm)
            np.save('Scratch/freqs', freqs)

        #  ##Write output files
        if self.write_output_files:
            outputFile = 'Output/{}_{}{}_{}.dat'.format(self.planet, self.outType, btmp, runStart.strftime("%Y%m%d_%H%M"))
            if self.verbose:
                print('\nWriting {} data to {}'.format(outType, datFile))
            self.__setHeader__(missed_planet)
            self.fIO.write(outputFile, outType, freqs, freqUnit, b, self.Tb, self.header)
            if self.plot and outType == 'Profile':
                plt.figure("Profile")
                Tbtr = np.transpose(self.Tb)
                for j in range(len(freqs)):
                    frqs = ('%.2f %s' % (self.freqs[j], self.freqUnit))
                    plt.plot(bs, Tbtr[j], label=frqs)
                plt.legend()
                plt.xlabel('b')
                plt.ylabel('$T_B$ [K]')

        return self.data_return

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

    def __setHeader__(self, missed_planet):
        if missed_planet:
            self.header['res'] = '# res not set\n'
            self.header['orientation'] = '# orientation not set\n'
            self.header['aspect'] = '# aspect tip, rotate not set\n'
            self.header['rNorm'] = '# rNorm not set\n'
        else:
            self.header['orientation'] = '# orientation:   {}\n'.format(repr(self.config.orientation))
            self.header['aspect'] = '# aspect tip, rotate:  {:.4f}  {:.4f}\n'.format(utils.r2d(self.tip), utils.r2d(self.rotate))
            self.header['rNorm'] = '# rNorm: {}\n'.format(self.rNorm)
            if self.bType:
                self.header['bType'] = '# bType:  {}\n'.format(self.bType)
            if self.outType:
                self.header['outType'] = '# outType:  {}\n'.format(self.outType)
                if self.outType == 'Image':
                    self.header['imgSize'] = '# imgSize: {:.0f}, {:.0f}\n'.format(self.imRow, self.imCol)
                    resolution = utils.r2asec(math.atan(abs(self.b[1][0] - self.b[0][0]) * self.rNorm / self.config.distance))
                    print('resolution = ', resolution)
                    self.header['res'] = '# res:  {} arcsec\n'.format(resolution)
        self.header['gtype'] = '# gtype: {}\n'.format(self.config.gtype)
        self.header['radii'] = '# radii:  {:.1f}  {:.1f}  km\n'.format(self.config.Req, self.config.Rpol)
        self.header['distance'] = '# distance:  {} km\n'.format(self.config.distance)

    def set_b(self, b, block):
        """b has a number of options for different bType:
               'points':  discrete number of points
               'line':  radial lines
               'image':  full image
               'stamp':  small image of region
               'disc':  disc-averaged
           outType can be 'image', 'spectrum', 'profile'
           b, bType and outType get set
           bType is never used anywhere..."""
        self.header['b'] = '# b request:  {}  {}\n'.format(str(b), str(block))

        self.imSize = None
        if type(b) == float:  # this generates a grid at that spacing and blocking
            bType = 'image'
            outType = 'Image'
            pb = []
            grid = -1.0 * np.flipud(np.arange(b, 1.5 + b, b))
            grid = np.concatenate((grid, np.arange(0.0, 1.5 + b, b)))
            # get blocks
            bsplit = len(grid) / abs(block[1])
            lastRow = block[0] / abs(block[1])
            if abs(block[1]) == 1:
                lastRow = 0
            for i in range(int(bsplit + lastRow)):
                ii = i + int((block[0] - 1) * bsplit)
                vrow = grid[ii]
                for vcol in grid:
                    pb.append([vcol, vrow])
            b = pb
            self.imSize = [len(grid), len(b) / len(grid)]
        elif type(b) == str:
            bType = b.lower()
            if bType not in ['disc', 'stamp']:
                raise ValueError('Invalid b request string')
            if bType == 'disc':
                b = [[0.0, 0.0]]
                outType = 'Spectrum'
                print('Setting to disc-average')
            elif bType == 'stamp':
                outType = 'Image'
                print('Setting to postage stamp.  Need more information')
                try:
                    bres = float(raw_input('...Input postage stamp resolution in b-units:  '))
                except ValueError:
                    bres = None
                bxmin, bxmax = raw_input('...Input bx_min, bx_max:  ').split(',')
                try:
                    bxmin = float(bxmin)
                    bxmax = float(bxmax)
                except ValueError:
                    bres = None
                bymin, bymax = raw_input('...Input by_min, by_max:  ').split(',')
                try:
                    bymin = float(bymin)
                    bymax = float(bymax)
                except ValueError:
                    bres = None
                if bres:
                    pb = []
                    for x in np.arange(bxmin, bxmax + bres / 2.0, bres):
                        for y in np.arange(bymin, bymax + bres / 2.0, bres):
                            pb.append([y, x])
                    b = pb
                    xbr = len(np.arange(bymin, bymax + bres / 2.0, bres))
                    self.imSize = [xbr, len(b) / xbr]
        elif len(np.shape(b)) == 1:     # this makes:
            pb = []
            if len(b) == 2:
                bType = 'points'
                outType = 'Spectrum'
                pb.append(b)            # ...data at one location
            elif len(b) > 3:
                bType = 'line'
                outType = 'Profile'
                outType = 'Spectrum'
                angle = utils.d2r(b[0])
                del b[0]
                if len(b) == 3:
                    b = np.arange(b[0], b[1] + b[2] / 2.0, b[2])  # ...this generates the points as start,stop,step
                for v in b:
                    pb.append([v * math.cos(angle), v * math.sin(angle)])  # ...a line at given angle (angle is first entry)
            else:
                raise ValueError('Invalid b request.')
            b = pb
        elif len(np.shape(b)) == 2:
            bType = 'points'
            if len(b) > 10:
                outType = 'Profile'
                outType = 'Spectrum'
            else:
                outType = 'Spectrum'
        else:
            raise ValueError('Invalid format for b.')

        if self.batch_mode:
            outType = 'Spectrum'
        # ##Set as self for header information
        self.outType = outType
        self.bType = bType

        return b, outType

    def set_freq(self, freqs, freqUnit):
        """ Internal processing of frequency list.
               if freqs is a string, it reads frequencies from that file
               if freqs is a scalar, it is made to a list of length 1
               if freqs is a list with first entry a string:
                    'r'ange:  sets range with 3 following values in list as start,stop,step
                        if the step is negative, it is assumed as a log step
                    that is currently the only option
               otherwise the list is used"""
        self.header['freqs'] = '# freqs request: {} {}\n'.format(str(freqs), freqUnit)
        # ## Process frequency range "request"
        if type(freqs) == str:
            freqs = np.loadtxt(freqs)
            freqs = list(freqs)  # convert numpy array back to list (for no good reason really)
        elif type(freqs) != list:
            freqs = [freqs]
        # We now have a list, check if first entry is a keyword
        if type(freqs[0]) == str:
            if freqs[0][0].lower() == 'r' and len(freqs) == 4:
                fstart = freqs[1]
                fstop = freqs[2]
                fstep = freqs[3]
                freqs = []
                f = fstart
                if (fstep < 0.0):  # do log steps
                    fstep = abs(fstep)
                    while f <= fstop:
                        freqs.append(f)
                        f *= fstep
                else:
                    while f <= fstop:
                        freqs.append(f)
                        f += fstep
            else:
                raise ValueError('Invalid format for frequency request')
        for i in range(len(freqs)):
            freqs[i] = utils.convert_unit(freqs[i], freqUnit)
        if len(freqs) > 1:
            s = '{} in {} frequency steps ({} - {} {})'.format(self.planet, len(freqs), freqs[0], freqs[-1], utils.proc_unit(freqUnit))
        else:
            s = '{} at {} {}'.format(self.planet, freqs[0], utils.proc_unit(freqUnit))
        if not self.super_quiet:
            utils.log(self.log, s, True)
        self.freqs = freqs
        self.freqUnit = utils.proc_unit(freqUnit)
        return freqs, utils.proc_unit(freqUnit)
