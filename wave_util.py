#!/bin/python2
import numpy as np
import math
import random

"""
zerodown takes elevation data and its sampling ratio and computes values for individual wave height and period sorted descending by wave height. In addition it returns values for significant(1/3) and 1/10 individual wave height/period.
"""
def zerodown(elevation, sampling_ratio):
    # make no assumptions regarding the unit of time
    time_length = len(elevation) * sampling_ratio

    # 
    iz = 0
    i_prev = 0

    #
    tz = 0.0
    tz_prev = 0.0

    ht = []

    # convert into numpy array
    elevation = np.asarray(elevation, float)
    # calculate mean sea level of data
    msl = -np.mean(elevation)
    # normalize elevation
    elevation = np.add(elevation, msl)

    # determine height and period of each individual wave
    for i in xrange(1,len(elevation)-1):
        # if elevation at current index is positive and elevation at the next index is negative
        if elevation[i] >= 0 and elevation[i+1] < 0:
            # calculate the interpolated time of zero crossing
            tz = lagrpol(elevation[i:i+2], xrange(i,i+2), 2, 0.0)
            # calculate period of current individual wave
            T = sampling_ratio * (tz - tz_prev)

            # slice elevation data to current individual wave
            Z = elevation[i_prev:i]
            # calculate height of current individual wave
            H = np.max(Z) - np.min(Z)

            # append to our ouput array
            ht.append([H, T])

            iz += 1
            tz_prev = tz
            i_prev = i

    # sort descending and convert to numpy array
    ht.sort()
    ht = ht[::-1]
    htsort = np.asarray(ht,float)

    # obtain mean, hs,ts and h10,t10 through running average pass
    n = len(htsort)
    ns = n/3
    n10 = n/10
    havg = 0
    tavg = 0
    for i in xrange(n):
        havg += (htsort[i,0] - havg) / (i+1)
        tavg += (htsort[i,1] - tavg) / (i+1)
        if i == n10:
            ht10 = [havg, tavg]
        if i == ns:
            hts = [havg, tavg]

    mean = [havg, tavg]

    hmax = htsort[0,0]
    tmax = np.max(htsort[:,1])
    max = [hmax, tmax]

    return (mean, max, hts, ht10)

"""
lagrange interpolation
"""
def lagrpol(x, y, n, xx):
    fx = 0.0
    l = 1.0

    for i in xrange(n):
        l = 1.0
        for j in xrange(n):
            if j != i:
                l *= (xx - x[j]) / (x[i] - x[j])
        fx += l * y[i]

    return fx

"""
generate jonswap spectrum
"""
def spectrum(hs, ts, sr, f):
    sp = []                      # array for our output spectrum

    # spectrum generator constants
    gamma_1 = 3.3
    beta_i = 0.0624*(1.094 - 0.01915*math.log10(gamma_1))/(0.230 + 0.0336*gamma_1 - (0.185/(1.9 + gamma_1)))
    divtstp = 1 - 0.132*(gamma_1 + 0.2)**(-0.559)
    divltfp = 2*0.07**2
    divgtfp = 2*0.09**2
    tp = ts/divtstp
    hssq = hs**2
    tppow4 = tp**(-4)

    fp = 1/tp
    for i in f:
        # two condition of sepctra depending on the value of fp
        if i <= fp:
            sp.append(beta_i*hssq*tppow4*i**(-5)*math.exp(-1.25*(i*tp)**(-4))*gamma_1**(math.exp(-(i*tp - 1)**2)/divltfp))
        else:
            sp.append(beta_i*hssq*tppow4*i**(-5)*math.exp(-1.25*(i*tp)**(-4))*gamma_1**(math.exp(-(i*tp - 1)**2)/divgtfp))

    return sp

"""
sea water elevation generation from spectrum
"""
def ema(sp, sr, f, t, del_f):
    eta = [0.0]*7200    # array containing our ouput
    rand = random.seed()
    pi2 = 2*math.pi
    a = np.asarray(S,float)
    a *= 2*del_f**0.5
    for i in np.arange(len(t)):
        eta[i] = sum((a[j]*math.cos(pi2*f[j]*t[i]+(random.random()*pi2))) for j in np.arange(len(f)))
    
