from __future__ import print_function, absolute_import, division
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os
import numpy as np

commentChars = ['!', '#', '$', '%', '&', '*']
affirmative = [1, '1', 'y', 'Y', 't', 'T']

Units = {'Frequency': 1, 'Hz': 1.0, 'kHz': 1.0E3, 'MHz': 1.0E6, 'GHz': 1.0E9,
         'Length': 1, 'm': 1.0, 'km': 1.0E3, 'cm': 1.0E-2, 'AU': 149597870691.0,
         'Pressure': 1, 'bars': 1.0, 'atm': 1.01325,
         'Time': 1, 'sec': 1.0, 'min': 60.0, 'hr': 3600.0, 'day': 86400.0, 'year': 31536000.0,
         'Acceleration': 1, 'mpersec2': 1.0, 'cmpersec2': 0.01}
processingUnits = {'GHz': ['GHz', 'Hz', 'kHz', 'MHz'],
                   'km': ['m', 'cm', 'AU', 'km'],
                   'bars': ['bars', 'atm'],
                   'sec': ['sec', 'min', 'hr', 'day', 'year'],
                   'mpersec2': ['mpersec2', 'cmpersec2']}


def proc_unit(supplied_unit):
    proc_unit = None
    for proc, units in processingUnits.iteritems():
        if supplied_unit in units:
            proc_unit = proc
            break
    return proc_unit


def convert_unit(v, supplied_unit):
    converted = v
    for proc, units in processingUnits.iteritems():
        if supplied_unit in units:
            converted = v * Units[supplied_unit] / Units[proc]
            break
    return converted


alphaUnit = 'invcm'
atmLayerUnit = 'km'
processingFreqUnit = proc_unit('Hz')
speedOfLight = 2.9979E8     	# m/s
kB = 1.3806503E-23          	# m2 kg s-2 K-1  Boltzmann's constant
hP = 6.626068E-34       	# m2 kg / s	 Planck's constant
T_cmb = 2.725

rfBands = {'HF': [0.003, 0.03], 'VHF': [0.03, 0.3], 'UHF': [0.3, 1.0], 'L': [1.0, 2.0], 'S': [2.0, 4.0],
           'C': [4.0, 8.0], 'X': [8.0, 12.0], 'Ku': [12.0, 18.0], 'K': [18.0, 26.5], 'Ka': [26.5, 40.0],
           'Q': [40.0, 50.0], 'V': [50.0, 75.0], 'W': [75.0, 110.0]}


def r2d(a):
    return a * 180.0 / np.pi


def d2r(a):
    return a * np.pi / 180.0


def r2asec(a):
    return 3600.0 * r2d(a)


def getRFband(freq, unit='GHz'):
    freq = freq * Units[unit] / Units['GHz']
    for bnd in rfBands:
        if rfBands[bnd][0] <= freq < rfBands[bnd][1]:
            return bnd
    return None


def invertDictionary(dic):
    e = {}
    for d in dic.keys():
        e[dic[d]] = d
    sk = e.keys()
    sk.sort()
    return e, sk


def setupLogFile(log):
    if type(log) == file:
        logfp = log
    elif type(log) == str:
        logfp = open(log, 'a')
    else:
        logfp = None
    return logfp


def log(logfp, msg, printOut=True):
    if logfp:
        logfp.write(msg + '\n')
        if printOut:
            print(msg)


def close(logfp):
    if logfp:
        logfp.close()


def ls(directory='Output', tag='dat', show=True, returnList=False):
    """Generates file list for plotTB and writeWavel"""
    filelist = os.listdir(directory)
    files = []
    i = 0
    for fff in filelist:
        show_tag = tag is None or (isinstance(tag, str) and tag in fff)
        if fff[0] != '.' and show_tag:
            files.append(os.path.join(directory, fff))
            if show:
                print('{}:  {}'.format(i, files[i]))
            i += 1
    if returnList:
        return files


def get_data_from(line):
    if line[0] in commentChars or len(line) < 4:
        return None
    data = line.split()
    try:
        cval = [float(x) for x in data]
    except ValueError:
        cval = None
    return cval


def get_expected_number_of_entries(fp):
    """Ad hoc function to guess the number of file entries expected
    """
    enoe = {}
    for line in fp:
        cval = get_data_from(line)
        if cval is not None:
            enoe[len(cval)] = enoe.setdefault(len(cval), 0) + 1
    fp.seek(0)
    vm = [-1, 0]
    for v in enoe.iteritems():
        if v[1] > vm[1]:
            vm = v
    return v[0]


def bang():
    files = ls(show=False, returnList=True)
    for f in files:
        plotTB(f, xaxis='wavel', xlog=True, justFreq=True)


def writeWavel(fn=None, outputFile=None, directory='Output'):
    filename, Tb, f, wavel, b, xlabel, ylabels = readTB(fn=fn, directory=directory)
    title = filename.split('/')[-1].split('.')[0]

    # Write file
    title += '_wavel.dat'
    if outputFile is None:
        outputFile = title
    print('Writing to ', outputFile)
    fp = open(outputFile, 'w')
    for i in range(len(wavel)):
        s = '%f\t' % (wavel[i])
        fp.write(s)
        for j in range(len(b)):
            s = '%f\t' % (Tb[j][i])
            fp.write(s)
        fp.write('\n')
    fp.close()
    return n


def concatdat(files, directory='Output'):
    """Given a list of utils.ls indices, returns TB, freq, wavel"""
    aTB = []
    bf = []
    cwavel = []
    for fil in files:
        data = readTB(fn=fil)
        a = list(data[1][0])
        b = list(data[2])
        c = list(data[3])
        aTB += a
        bf += b
        cwavel += c
    sa = np.argsort(np.array(bf))
    TB = []
    f = []
    wavel = []
    for i in sa:
        TB.append(aTB[i])
        f.append(bf[i])
        wavel.append(cwavel[i])
    return TB, f, wavel


def plotObs(fn, cols=[0, 1, 2], color='b', marker='o', delimiter=None, comline='!'):
    try:
        fp = open(fn, 'r')
    except IOError:
        print(fn, ' not found')
        return 0
    data = []
    for line in fp:
        if comline in line[0:5]:
            continue
        dline = line.split(delimiter)
        if len(dline) < max(cols):
            continue
        drow = []
        for c in cols:
            drow.append(float(dline[c]))
        data.append(drow)
    data = np.array(data)
    plt.semilogx(data[:, 0], data[:, 1], color=color, marker=marker)
    plt.errorbar(data[:, 0], data[:, 1], yerr=data[:, 2], color=color, marker=marker, ls='none')
