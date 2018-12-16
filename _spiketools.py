import math
import random

import numpy as np


def poissonTrain(rate, tmax, tstart=4., seed=None):
    '''rate [kHz], tmax [ms], tstart [ms]'''
    if rate > 1e-12:
        if seed != None:
            random.seed(seed)
        spiketrain = []
        p = -math.log(1.0 - random.random()) / rate + tstart
        while p < tmax:
            spiketrain.append(p)
            p = (-math.log(1.0 - random.random()) / rate) + spiketrain[-1]
        return spiketrain
    else:
        return []


def modulatedPoissonTrain(rate, modrate, tmax, tstart=2.):
    '''rate [kHz], modrate [kHz], tmax [ms], tstart [ms]'''
    pass


def piecewisePoissonTrain(rates=[], ts=[], tstart=1., seed=None):
    ''' rates [kHz], ts [ms], tstart[ms] '''
    if seed != None:
        random.seed(seed)
    spiketrain = []
    t0 = tstart
    for i, rate in enumerate(rates):
        t1 = ts[i]
        spktr = poissonTrain(rate, t1, t0)
        spiketrain.extend(spktr)
        t0 = t1

    return spiketrain


def inhomogeneousPoissonTrain(rate_arr, dt, tstart=1.):
    indicator_arr = np.random.rand(*rate_arr.shape)
    return dt * np.where(indicator_arr < rate_arr*dt)[0] + tstart


def spikesFromVm(t, vm, **kwargs) :
    """ Extract the spike times from a voltage trace.

    Parameters
    ----------
    t : array or list of int or float
        time points (in ms)
    vm : array or list of int or float
        voltage points

    threshold : int or float, optional
        voltage at which we consider a spike to be fired
    min_t : int or float, optional
        only spikes after this moment are considered
    max_t : int or float, optional
        only spikes before this moment are considered

    Returns
    -------
    spiketimes : list of int or float
        spike times (in ms)
    """
    T = kwargs['threshold'] if('threshold' in kwargs) else 0
    min_t = kwargs['min_t'] if('min_t' in kwargs) else 0
    max_t = kwargs['max_t'] if('max_t' in kwargs) else max(t)
    in_spike = False
    spiketimes = []
    for i in range(len(t)):
        if (vm[i] > T) and (not in_spike): # SLOW
            if min_t <= t[i] <= max_t :
                spiketimes.append(t[i])
            in_spike = True
        elif( (in_spike == True) and (vm[i] < T) ):
            in_spike = False
    '''
    # if you want, do a visual check with matplotlib
    plt.plot(data[:,0],data[:,1])
    plt.plot(spiketimes[:],np.tile(0,len(spiketimes) ), 'rv'  )
    plt.savefig('temp_spikesplot.pdf')
    plt.show()
    '''
    return spiketimes