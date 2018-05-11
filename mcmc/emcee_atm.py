import sys
import os
import emcee
import corner
from shutil import copyfile
import numpy as np
import time
import mcmc.emcee_input as emcee_input
import planet
import matplotlib.pyplot as plt

gen = None
# ## Load in real data and model spectrum


def get_obs_spectrum(specdatfile):
    '''Read spectrum data file
        from ./Output/.dat
        Should have format:
        wavelengths (GHz), Tb (K), Tb_err(K)'''

    spec = []
    with open(specdatfile, 'r') as tmp:
        for line in tmp:
            if not line.startswith("#"):
                values = [float(x) for x in line.split()]
                if len(spec) == 0:
                    spec = values
                else:
                    spec = np.vstack((spec, values))

    return spec


def get_model_spectrum(p, val, freqs, b, par_names, limits):
    '''Generate model spectrum'''

    new_scale_factors = update_scale(p, par_names, val, limits)
    for c in par_names:
        p.alpha.scale_constituent_values[c] = new_scale_factors[c]

    data = p.run(freqs=freqs, b=b)

    spec = np.vstack((data.f, data.Tb))
    spec = spec.T

    return spec


def update_scale(p, par_names, guess, limits):
    '''Update scale val with new val'''

    return gen.generate(p, par_names, guess)

# ## Compute probabilities


def lnprior(theta, limits):  # flat priors
    tmp = 1
    for i in range(len(theta)):
        if limits[i][0] < theta[i] < limits[i][1]:
            tmp *= 1
        else:
            tmp *= 0
    if tmp == 1:
        return 0.0
    else:
        return -np.inf


def lnlike(theta, x, y, yerr, p, freqs, b, par_names, limits):
    parvals = theta
    spec = get_model_spectrum(p, parvals, freqs, b, par_names, limits)
    ymodel = spec[:, 1]

    sigsq = yerr**2
    lnP = -0.5 * np.sum((y - ymodel)**2. / sigsq + np.log(2 * np.pi * sigsq))
    return lnP


def lnprob(theta, x, y, yerr, p, freqs, b, par_names, limits):
    lp = lnprior(theta, limits)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr, p, freqs, b, par_names, limits)


def run_emcee_spectrum(sampler, pos, nsteps, outdatfile, lnprobfile=None):
    t0 = time.time()
    counter = 0
    print('Starting MCMC')
    for result in sampler.sample(pos, iterations=nsteps, storechain=True):
        counter += 1
        if counter % 10 == 0:
            print('Iteration %d' % counter)
        position = result[0]
        f = open(outdatfile, "a")
        # pdb.set_trace()
        for k in range(position.shape[0]):
            f.write("{0:4d} {1:s}\n".format(k, " ".join(map(str, position[k]))))
        f.close()
        if lnprobfile is not None:
            probability = result[1]
            f = open(lnprobfile, "a")
            for k in range(probability.shape[0]):
                f.write("{0:4d} {1:s}\n".format(k, repr(probability[k])))
            f.close()
    t1 = time.time()
    print ('run took {0} minutes to travel {1} steps'.format((t1 - t0) / 60., nsteps))
    print ('acceptance fraction:', sampler.acceptance_fraction)

    return sampler


def run_emcee_spectrum_new(config='config.par', output=True):
    global gen
    inp = emcee_input.gen_emcee_input()

    name = inp['planet']
    refData = inp['refData']

    freqs = inp['freqs']
    b = inp['b']

    par_names = inp['parameters']['names']
    guess = inp['parameters']['guesses']
    limits = inp['parameters']['limits']
    nwalker = inp['nwalkers']
    threads = inp['threads']
    nsteps = inp['nsteps']

    outdatfile = inp['outdatafile']
    lnprobfile = inp['lnprobfile']
    configfile = inp['configfile']

    # generate the absorb.npy files once and import scale compute module
    p = planet.Planet(name, config=configfile, generate_alpha='True', plot=False)
    p.run(freqs=freqs, b=b)
    sys.path.append(p.planet)
    __import__(p.config.scalemodule)
    gen = sys.modules[p.config.scalemodule]

    # initialize scale file if needed
    if not os.path.isfile(p.config.scale_file_name):
        print("{} does not existing - initializing new one".format(p.config.scale_file_name))
        atm_pressure = p.atm.gas[p.atm.config.C['P']]
        planet.alpha.initialize_scalefile(p.config.scale_file_name, atm_pressure, par_names, guess)

    # setup for mcmc
    p = planet.Planet(name, mode='mcmc', config=configfile)

    # Check if we're going to overwrite a file
    if os.path.isfile(outdatfile):
        overwrite = True
        if overwrite is False:
            print ("Please choose another filename and rerun: {}".format(outdatfile))
            return
    continueprob, newprob = True, True  # = check_likefile(lnprobfile, "lnprobfile")
    if continueprob is False:
        print ("Exiting because of a problem with lnprobfile.")
        return
    if newprob is False:
        print ("lnprobfile already exists! Exiting")
        return
    f = open(lnprobfile, "w")
    f.close()

    name = str.capitalize(name)
    ndim = len(guess)

    obs_spec = get_obs_spectrum(refData)
    x = obs_spec[:, 0]
    y = obs_spec[:, 1]
    yerr = obs_spec[:, 2]

    pos = [guess + 1e-2 * np.random.randn(ndim) for i in range(nwalker)]
    sampler = emcee.EnsembleSampler(nwalker, ndim, lnprob, args=(x, y, yerr, p,
                                    freqs, b, par_names, limits), threads=threads)

    f = open(outdatfile, "w")
    f.close()
    sampler = run_emcee_spectrum(sampler, pos, nsteps, outdatfile, lnprobfile=lnprobfile)

    # Write out and plot scale
    if output:
        output_scale = 'Scratch/scale_mcmc_out.dat'
        p.alpha.write_scale(output_scale)
        p.atm.scaleAtm(output_scale, plot_diff=True)

    return sampler


def run_emcee_spectrum_append(output=True):
    global gen
    inp = emcee_input.gen_emcee_input()

    name = inp['planet']
    refData = inp['refData']

    par_names = inp['parameters']['names']
    guess = inp['parameters']['guesses']
    limits = inp['parameters']['limits']
    nwalker = inp['nwalkers']
    threads = inp['threads']
    nsteps = inp['nsteps']

    freqs = inp['freqs']
    b = inp['b']

    datfile = inp['outdatafile']
    lnprobfile = inp['lnprobfile']
    configfile = inp['configfile']
    initial_scalefile = inp['scalefile']

    name = name.capitalize()

    obs_spec = get_obs_spectrum(refData)
    x = obs_spec[:, 0]
    y = obs_spec[:, 1]
    yerr = obs_spec[:, 2]

    if os.path.isfile(initial_scalefile) is False:
        print ("{0} is not found. Did you mean to run run_emcee_spectrum_new?".format(initial_scalefile))
        return
    if os.path.isfile(datfile) is False:
        print ("{0} is not found. Did you mean to run run_emcee_spectrum_new?".format(datfile))
        return
    append = True
    print ("appending to file {0}?".format(datfile))
    if append is False:
        print ("If you want to continue running emcee on {0}, then say so! Returning.".format(datfile))
        return
    # Check probability file
    continueprob, newprob = check_likefile(lnprobfile, "lnprobfile")
    if continueprob is False:
        print ("Exiting because of a problem with lnprobfile.")
        return
    if newprob is True:
        print ("Warning: printing probabilities to a NEW file, while appending walker positions to an OLD file.")
        f = open(lnprobfile, "w")
        f.close()

    samples, ndim, nwalkers = read_emcee_datfile(datfile)
    pos = samples[-1, :, :]
    sampler = emcee.EnsembleSampler(nwalker, ndim, lnprob, args=(x, y, yerr, p, freqs, b, par_names, limits), threads=threads)
    sampler = run_emcee_spectrum(sampler, pos, nsteps, datfile, lnprobfile=lnprobfile)
    # Write out and plot scale
    if output:
        output_scale = 'Scratch/scale_mcmc_out.dat'
        p.alpha.write_scalefile(output_scale)
        p.atm.scaleAtm(output_scale, plot_diff=True)

    return sampler


def check_likefile(filename, filestring):
    new = False
    if filename is None:
        continue_run = False
    else:
        if os.path.isfile(filename) is False:
            print ('writing new file')
            continue_run = True
            new = True
        else:
            continue_run = True
            print ('appending to existing file')
    return continue_run, new


def read_emcee_datfile(datfile):
    data = np.loadtxt(datfile).T
    nwalkers = np.max(data[0]) + 1
    niter = data[0].shape[0] / nwalkers
    ndim = data.shape[0] - 1
    samples = data[1:].reshape((int(nwalkers), int(niter), int(ndim)))
    return samples


def read_emcee_probfile(probfile):
    with open(probfile, 'r') as file:
        tmp = file.readlines()
    out = []
    for line in tmp:
        values = [float(x) for x in line.split()]
        if len(out) == 0:
            out = values
        else:
            out = np.vstack((out, values))

    nwalkers = np.max(out[:, 0]) + 1
    probabilities = out[:, 1].reshape(-1, nwalkers)
    return probabilities, nwalkers


def plot_burnin(samples, names, cut=0):
    nwalkers = samples.shape[0]
    niter = samples.shape[1]
    ndim = samples.shape[2]

    fig, axes = plt.subplots(1, ndim, figsize=(5, 5 * ndim))

    for i in range(ndim):
        label = names[i]
        if ndim == 1:
            ax = axes
        else:
            ax = axes[i]

        for w in range(nwalkers):
            ax.plot(np.arange(1, niter + 1), samples[w, :, i])
        ax.set_ylabel(label)
        ax.set_xlabel('Iteration')

        flatchain = samples[:, cut:, i].flatten()
        [lower, middle, upper] = np.percentile(flatchain, [16, 50, 84])
        ax.axhline(middle, color='k', linestyle='--')
        ax.axhline(lower, color='k', linestyle=':')
        ax.axhline(upper, color='k', linestyle=':')

    plt.tight_layout()
    plt.show()


def best_fit(samples, labels, cut=0):
    '''Given a chain of MCMC samples, output 16th, 50th, and 84th percentile
    in each dimension. labels has same length as ndim'''

    ls = []
    ms = []
    us = []
    for i in range(samples.shape[2]):
        flatchain = samples[:, cut:, i].flatten()
        [lower, middle, upper] = np.percentile(flatchain, [16, 50, 84])
        ls.append(lower)
        ms.append(middle)
        us.append(upper)
    print('Labels', labels)
    print('Best Fits', ms)
    print('Lower Error', ls)
    print('Upper Error', us)
    return ms, ls, us


def cornerplot(samples, labels, cut=0):
    samples = samples[:, cut:, :]
    fig = corner.corner(samples.reshape((-1, samples.shape[2])), labels=labels)
    # truths = [x1_true, x2_true, ...)

    plt.show()

# def best_fit(samples, names, burnin):
#     flatchain=samples.swapaxes(0,1)[:,burnin:,:].reshape(-1,samples.shape[2])
#     results  = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
#                    zip(*np.percentile(flatchain, [16, 50, 84], axis=0)))
#     for i in range(len(names)):
#         print '{0} : {1} +{2}, -{3}'.format(names[i],results[i][0],results[i][1],results[i][2])
#     return results
#
# def plot_simple_burnin(samples, nwalkers):
#     print(samples.shape)
#     for i in range(nwalkers):
#         plt.plot(samples[:,i,0])
#
#     plt.show()
