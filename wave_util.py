#!/bin/python2
import numpy as np
import numpy.random as nprand
import numpy.fft as npf
import numexpr as ne
import math

"""
zerodown takes elevation data and its sampling ratio and computes values for individual wave height and period sorted descending by wave height. In addition it returns values for significant(1/3) and 1/10 individual wave height/period.
"""
def zerodown(elevation, sampling_ratio):
    # make no assumptions regarding the unit of time
    srd = 1.0/sampling_ratio
    time_length = len(elevation) * srd

    # 
    i_prev = 0

    #
    t = 0.0
    t_prev = 0.0

    ht = []

    # lagrange polynomial dimension
    n = np.arange(2)

    # normalize elevation
    elevation -= np.mean(elevation)

    # determine height and period of each individual wave
    for i in xrange(1,len(elevation)-1):
        # if elevation at current index is positive and elevation at the next index is negative
        if elevation.item(i) >= 0 and elevation.item(i+1) < 0:
            # calculate the interpolated time of zero crossing
            t = lagrpol(elevation[i:i+2], np.arange(i,i+2), n, 0.0)
            # calculate period of current individual wave
            T = srd * (t - t_prev)

            # calculate height of current individual wave
            H = np.ptp(elevation[i_prev:i])

            # append to our ouput array
            ht.append([H, T])

            t_prev = t
            i_prev = i

    # sort descending and convert to numpy array
    ht.sort()
    ht[:] = ht[::-1]
    htsort = np.asarray(ht,float)

    # obtain mean, hs,ts and h10,t10 through running average pass
    n = len(htsort)
    n10 = n/10
    ns = n/3

    ht10 = [np.mean(htsort[:n10,0]), np.mean(htsort[n10:,1])]
    hts = [np.mean(htsort[:ns,0]), np.mean(htsort[:ns,1])]
    mean = [np.mean(htsort[:,0]), np.mean(htsort[:,1])]
    tmax = np.max(htsort[:,1])
    max = [htsort[0,0], tmax]

    """
    havg = 0
    tavg = 0
    for i in xrange(n):
        havg += (htsort.item(i,0) - havg) / (i+1)
        tavg += (htsort.item(i,1) - tavg) / (i+1)
        if i == n10:
            ht10 = [havg, tavg]
        if i == ns:
            hts = [havg, tavg]

    mean = [havg, tavg]

    tmax = np.max(htsort[:,1])
    max = [htsort.item(0,0), tmax]
    """

    return (mean, max, hts, ht10)


"""
lagrange interpolation
"""
def lagrpol(x, y, n, xx):
    fx = 0.0

    for i in n:
        l = 1.0
        for j in n:
            if j != i:
                l *= (xx - x.item(j)) / (x.item(i) - x.item(j))
        fx += l * y.item(i)

    return fx


"""
generate jonswap spectrum
"""
def spectrum(hs, ts, gamma_1, beta_i, divtstp, divltfp, divgtfp, max_t, f):
    sp = []     # array for our output spectrum

    # spectrum generator constants
    tp = ts/divtstp
    hssq = hs**2
    tppow4 = tp**(-4)

    # vanilla implementation
    """
    fp = 1.0/tp
    for i in f:
        # two condition of sepctra depending on the value of fp
        if i <= fp:
            sp.append(beta_i*hssq*tppow4*i**(-5)*math.exp(-1.25*(i*tp)**(-4))*gamma_1**(math.exp(-((i*tp - 1)**2)/divltfp)))
        else:
            sp.append(beta_i*hssq*tppow4*i**(-5)*math.exp(-1.25*(i*tp)**(-4))*gamma_1**(math.exp(-((i*tp - 1)**2)/divgtfp)))
    """

    # calculate fp (frekeunsi pemisah) and figure out which frequency index it's closest to
    fp = int(max_t/tp)
    # split freqeuncy array into two separate subarrays at fp index
    flt, fgt = np.split(f,[fp])

    # calculate spectrum for each frequency subarray and combine them back: for frequencies less than fp (flt) use divltfp, otherwise (fgt) use divgtfp. in addition, use numexpr to parallelize non-array portions of the expression
    sp.extend(ne.evaluate('beta_i*hssq*tppow4*flt**(-5)*exp(-1.25*(flt*tp)**(-4))*gamma_1**(exp(-((flt*tp - 1)**2)/divltfp))'))
    sp.extend(ne.evaluate('beta_i*hssq*tppow4*fgt**(-5)*exp(-1.25*(fgt*tp)**(-4))*gamma_1**(exp(-((fgt*tp - 1)**2)/divgtfp))'))

    return sp


"""
sea water elevation generation from spectrum
"""
def ema(sp, sr, del_f):
    pi2 = 2*math.pi
    srd = 1.0/sr

    # vectorize spectrum and weigh
    a = (np.asarray(sp,float)*2*del_f)**0.5

    # vanilla implementation: slow as fuck!
    """
    eta = []
    for i in t:
        temp = 0.0
        for j in xrange(len(f)):
            temp += a[j] * math.cos(pi2*f[j]*t[i] + pi2*nprand.random_sample())
        eta.append(temp)
    """

    # list comprehension: x2 improvement still slow!
    """
    eta = [sum([a.item(j)*math.cos(pi2*(f.item(j)*i*srd + nprand.random_sample())) for j in xrange(len(f))]) for i in t]
    """

    # with vectorization: ~60x improvement a lot better
    """
    #eta = np.sum((a*np.cos(pi2*(f*t*srd + nprand.random_sample(f.shape)))), 0)
    """

    # with FFT: ~6000x improvement (!!!)
    eta = npf.irfft(a * np.exp(1j * nprand.uniform(0, pi2, a.shape))) * a.size

    return eta
