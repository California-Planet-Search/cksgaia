import numpy as np
import pandas as pd
import scipy.stats
import scipy.ndimage
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline, UnivariateSpline, LSQUnivariateSpline
import pp

import lmfit

import cksgaia.misc
from cksgaia.config import *


# Period and radius bins
Redges = np.logspace(np.log10(0.5), np.log10(20), 36)
Pedges = np.logspace(np.log10(1.0), np.log10(100), 18)

# Error scaling factor from simulations
efudge = np.array([ 2.81570765,  2.49756761,  2.30114599,  2.5389933 ,  2.34539564,
        2.0856459 ,  1.92438544,  1.94995261,  1.88536633,  1.45970263,
        1.65253424,  1.8102754 ,  1.38489473,  1.50448697,  1.3923987 ,
        1.58201561,  1.47604433,  1.57615124,  1.25218412,  1.48151002,
        1.47226211,  1.45848035,  1.62620336,  1.45336104,  1.49927148,
        1.51589925,  1.33522687,  1.43910017,  1.45971532,  1.51720192,
        1.56983652,  1.36153125,  1.34642358,  1.44756085,  1.4414092 ]) * 0 + 1.0
efudge = 1.0

# node locations for spline fit
nodes= np.log10([1.3, 1.5, 1.9, 2.4, 3.0, 4.5, 11])
#nodes= np.log10([1.4, 1.8, 2.4, 3.0, 4.0])

# mask limit for fitting
masklim = (1.14, 20.0)

Rbinsize = np.diff(np.log10(Redges)).mean()
mask2 = np.ones_like(Redges, dtype=bool)
mask2[np.where((Redges <= masklim[0]))[0]] = False
#mask2[np.where((Redges >= masklim[1] + 2*Rbinsize) | (Redges <= masklim[0]))[0]] = False


def gamma_complete(x, result=None, step=False):
    if result is not None:
        a = result.params['a'].value
        loc = result.params['loc'].value
        scale = result.params['scale'].value
    else:
        a = 17.56
        loc = 1.00
        scale = 0.49

    f = scipy.stats.gamma.cdf(x, a, loc, scale)
        
    if step:
        cutoff = 12.0
        if len(x) < 1:
            return int(f > cutoff)
        else:
            pos = x <= cutoff
            f[pos] = 0.0
            f[~pos] = 1.0

    return f


def logistic(x, a=0.0, b=8.06, c=8.11, d=0.995, step=False):
    """
    Christiansen et al. (2016) logistic
    function for Kepler completeness as a function of 
    MES.

    Args:
        x (float): MES at which to calculate completeness
    Returns
        float: completeness at MES=x
    """
    
    # a = 0.0
    # b = 8.06
    # c = 8.11
    # d = 0.995
    f = ((a - d)/(1+(x/c)**b)) + d

    if step:
        cutoff = 12.0
        if len(x) < 1:
            return int(f > cutoff)
        else:
            pos = x <= cutoff
            f[pos] = 0.0
            f[~pos] = 1.0

    return f

def lognorm(x, amp, mu, sig, unc_limit=0.0967):
    #pdf = amp/(x*sig*np.sqrt(2*np.pi)) * np.exp(-(((np.log(x) - np.log(mu))**2)/(2*sig**2)))
    #pdf = amp * np.log10(np.e)/(x*sig*np.sqrt(2*np.pi)) * np.exp(-(((np.log10(x) - np.log10(mu))**2)/(2*sig**2)))
    
    wid = np.sqrt(sig**2 + unc_limit**2) / np.log(10.)
    pdf = gauss(np.log10(x), amp, np.log10(mu), wid)
    return pdf


def twolognorm(x, amp1, mu1, sig1, amp2, mu2, sig2, unc_limit=0.0967):
    pdf = lognorm(x, amp1, mu1, sig1, unc_limit=unc_limit) + lognorm(x, amp2, mu2, sig2, unc_limit=unc_limit)
    
    return pdf

def threefunc(x, c1, c2, tau, amp, mu, sig, b1=2.5, unc_limit=0.0967):
    mod = np.zeros_like(x)
    
    
    reg2 = np.where(x <= b1)[0]
    ng = lognorm(x[reg2], amp, mu, sig, unc_limit=0)
    
    mod[reg2] = c1 - ng
    
    break_val = c1 - lognorm(b1, amp, mu, sig, unc_limit=0)
        
    reg1 = np.where(x > b1)[0]
    mod[reg1] = lognorm(x[reg1], break_val-c2, b1, tau, unc_limit=0) + c2

    if unc_limit > 0:
        sigma = unc_limit * x.mean() / np.diff(x).mean()
        mod = scipy.ndimage.gaussian_filter(mod, sigma)

    mod = np.clip(mod, 0, 1)
    
    return mod

def splinefunc(x, n1, n2, n3, n4, n5, n6, n7, unc_limit=0.0967):
    heights = [n1, n2, n3, n4, n5, n6, n7]
    
    mod = InterpolatedUnivariateSpline(nodes, heights, ext=3, k=2)(np.log10(x))

    if unc_limit > 0:
        #sigma = unc_limit * x.mean() / np.diff(x).mean()
        mu = 10**np.mean(np.log10(x))
        kernel = lognorm(x, 1.0, mu, unc_limit)
        kernel /= np.trapz(kernel)

        mod = scipy.ndimage.filters.convolve1d(mod, kernel)
        #mod = scipy.ndimage.gaussian_filter(mod, sigma)

    mod = np.clip(mod, 0, 1)
        
    return mod

def bin_model(xmod, ymod, xbins):

    bin_mod = []
    for i,b in enumerate(xbins[:-1]):
        inbin = np.where((xmod >= b) & (xmod < xbins[i+1]))[0]
        #print b, xbins[i+1], ymod[inbin].mean()
        bin_mod.append(ymod[inbin].mean())

    bin_mod = np.array(bin_mod)
        
    return bin_mod


def resid_spline(params, x, y, err, xbins):
    n1 = params['n1'].value
    n2 = params['n2'].value
    n3 = params['n3'].value
    n4 = params['n4'].value
    n5 = params['n5'].value
    n6 = params['n6'].value
    n7 = params['n7'].value

    mod = splinefunc(x, n1, n2, n3, n4, n5, n6, n7)
    bin_mod = bin_model(x, mod, xbins)
    
    resid = (y-bin_mod)/err
    #resid = (y-mod)/err

    resid[np.isnan(resid)] = 0
        
    return resid


def resid_piece(params, x, y, err):
    c1 = params['c1'].value
    c2 = params['c2'].value
    tau = params['tau'].value
    amp = params['amp'].value
    mu = params['mu'].value
    sig = params['sig'].value
    b1 = params['b1'].value
    
    mod = threefunc(x, c1, c2, tau, amp, mu, sig, b1)

    return ((y-mod)/err)


def resid_single(params, x, y, err):
    amp = params['amp'].value
    mu = params['mu'].value
    sig = params['sig'].value
    
    mod = lognorm(x, amp, mu, sig, unc_limit=0)
    res = ((y - mod)/err)
        
    #res = -0.5 * np.exp(np.sum(res**2))
        
    return res

def resid_two(params, x, y, err):
    amp1 = params['amp1'].value
    mu1 = params['mu1'].value
    sig1 = params['sig1'].value
    amp2 = params['amp2'].value
    mu2 = params['mu2'].value
    sig2 = params['sig2'].value

    mod = twolognorm(x, amp1, mu1, sig1, amp2, mu2, sig2, unc_limit=0)
    
    res = ((y - mod)/err)
    
    #res = -0.5 * np.exp(np.sum(res**2))
    
    return res


def gauss(x, amp, mu, sig, normed=False):
    m = amp*np.exp(-(x-mu)**2/(2*sig**2))
    if np.isnan(m).any():
        print amp, mu, sig
    
    if normed:
        m = m/np.trapz(m, x)
    
    return m

def twogauss(x, amp1, mu1, sig1, amp2, mu2, sig2):
    return gauss(x, amp1, mu1, sig1) + gauss(x, amp2, mu2, sig2)

# def comp_unc(physmerge, nsamp=10):
#     pm = physmerge.copy()
#
#     kic = cksrad.io.load_kic(filtered=False)
#     kicselect = cksrad.io.load_kic()
#     kicselect = cksrad.completeness.fit_cdpp(kicselect)
#
#     bin_centers = 10**(np.log10(Redges[:-1]) + np.diff(np.log10(Redges)).mean()/2)
#     Pcen = 10**(np.log10(Pedges[:-1]) + np.diff(np.log10(Pedges)).mean()/2)
#
#     hist_stack = []
#     for i in range(nsamp):
#         prad_noise = np.random.normal(0, pm['giso_prad_err'].values, size=len(pm))
#         pm['giso_prad'] += prad_noise
#
#         pm = cksrad.completeness.get_weights(pm, kicselect)
#
#         wdetections,_,_ = np.histogram2d(pm['koi_period'], pm['giso_prad'], bins=[Pedges,Redges], weights=pm['weight'])
#
#         hist_stack.append(wdetections)
#
#     hist_stack = np.array(hist_stack).std(axis=0)
#     err_hist = np.sqrt(np.sum(hist_stack**2, axis=0)) / num_stars
#
#     return err_hist

def histfit(physmerge, verbose=True, completeness=True, boot_errors=False):
    xbins = Redges
    
    if completeness:

        bin_centers = 10**(np.log10(Redges[:-1]) + Rbinsize/2)
        Pcen = 10**(np.log10(Pedges[:-1]) + np.diff(np.log10(Pedges)).mean()/2)

        detections,_,_ = np.histogram2d(physmerge['koi_period'], physmerge['giso_prad'], bins=[Pedges,Redges])
        wdetections,_,_ = np.histogram2d(physmerge['koi_period'], physmerge['giso_prad'], bins=[Pedges,Redges], weights=physmerge['weight'])

        if boot_errors:
                err_hist = comp_unc(physmerge)
        
        sfac = 1/physmerge['tr_prob'].mean()
        rhist = np.sum(detections, axis=0) 
        rhistn = rhist / num_stars * sfac
        rerr = np.sqrt(rhist) / num_stars * sfac

        whist = np.sum(wdetections, axis=0) 
        whistn = whist / num_stars
        if boot_errors:
            werr = np.sqrt((rerr * (whistn/rhistn) * efudge)**2 + err_hist**2)
        else:
            werr = rerr * (whistn/rhistn) * efudge

        N = whistn
        e = werr
    else:
        N, bin_edges = np.histogram(physmerge['giso_prad'], bins=xbins)
        bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])
        e = np.sqrt(N)

    mask = np.ones_like(bin_centers, dtype=bool)
    mask[np.where((bin_centers <= masklim[0]))[0]] = False
    #mask[np.where((bin_centers <= masklim[0]) & (bin_centers >= masklim[1])[0]] = False
    
    e[e==0] = np.inf
    e[~mask] = np.inf
    e[np.isnan(e)] = 100

    x = bin_centers[mask]
    y = N[mask]
    err = e[mask]

    err[0] = 10

    params = lmfit.Parameters()
    params.add('amp', value = 0.08, min=0.0, max=1.0)
    params.add('mu', value=2.0, min=0.01, max=10.0)
    params.add('sig', value=0.5, min=0.01, max=2.0)

    result = lmfit.minimize(resid_single, params, args=(x,y,err))
    if verbose:
        print "\n***************************"
        print "Single Gaussian fit"
        print lmfit.fit_report(result)

    params = lmfit.Parameters()
    params.add('amp1', value = 0.07, min=0.0, max=1.0, vary=True)
    params.add('mu1', value=1.29, min=0.01, max=10.0, vary=True)
    params.add('sig1', value=0.05, min=0.01, max=2.0, vary=True)
    params.add('amp2', value = 0.07, min=0.0, max=1.0, vary=True)
    params.add('mu2', value=2.38, min=0.01, max=10.0, vary=True)
    params.add('sig2', value=0.10, min=0.01, max=2.0, vary=True)

    result2 = lmfit.minimize(resid_two, params, args=(x,y,err))
    if verbose:
        print "\n***************************"
        print "Double Gaussian fit"
        print lmfit.fit_report(result2)

        
    params = lmfit.Parameters()
    params.add('c1', value = 0.07, vary=True)
    params.add('c2', value=0.01, vary=True)
    params.add('tau', value=0.10, min=0.01, max=2.0, vary=True)
    params.add('amp', value = 0.07, min=0.01, max=1.0, vary=True)
    params.add('mu', value=1.75, min=1.1, max=3.0, vary=True)
    params.add('sig', value=0.07, min=0.01, max=1.0, vary=True)
    params.add('b1', value=2.4, min=2.0, max=3.5, vary=True)
    
    result = lmfit.minimize(resid_piece, params, args=(x,y,err))
    if verbose:
        print "\n***************************"
        print "Piecewise fit"
        print lmfit.fit_report(result)
    

    params = lmfit.Parameters()
    params.add('n1', value = y[0], min=0.0, max=1.0, vary=False)
    params.add('n2', value=0.03, min=0.0, max=1.0, vary=True)
    params.add('n3', value=0.05, min=0.0, max=1.0, vary=True)
    params.add('n4', value = 0.10, min=0.0, max=1.0, vary=True)
    params.add('n5', value = 0.04, min=0.0, max=1.0, vary=True)
    params.add('n6', value = 0.004, min=0.0, max=1.0, vary=True)
    params.add('n7', value = 0.0005, min=0.0, max=1.0, vary=True)
    
    result = lmfit.minimize(resid_spline, params, args=(x,y,err,Redges[mask2]))
    if verbose:
        print "\n***************************"
        print "Spline fit"
        print lmfit.fit_report(result)

        
    return (mask, bin_centers, N, e, result, result2)


def mcmc(x, y, err, result, nsteps=1000):
    def lnprob_single(params):
        res = resid_single(params, x, y, err)

        return -0.5*np.sum(res**2 + np.log(2*np.pi*err**2))

    def lnprob_two(params):
        res = resid_two(params, x, y, err)
        
        return -0.5*np.sum(res**2 + np.log(2*np.pi*err**2))

    def lnprob_spline(params):
        res = resid_spline(params, x, y, err, Redges[mask2])

        ln = -0.5*np.sum(res**2 + np.log(2*np.pi*err**2))
        
        return ln
    
    mini = lmfit.Minimizer(lnprob_spline, result.params)
    nburn = int(np.round(nsteps/10.))
    print "Running %d steps after %d burn-in steps." % (nsteps, nburn)
    res = mini.emcee(burn=nburn, steps=nsteps, params=result.params, nwalkers=20)

    return res


def kde(values, errors):

    minloc = np.argmin(values)
    minlim = values[minloc] - 3*errors[minloc]
    
    maxloc = np.argmax(values)
    maxlim = values[maxloc] + 3*errors[maxloc]
    
    kx = np.linspace(minlim, maxlim, 2000)
    #kx = np.logspace(np.log10(minlim), np.log10(maxlim), 2000)
    ky = np.zeros_like(kx)
        
    for val,err in zip(values, errors):
        mod = gauss(kx, 1.0, val, err, normed=True)            
        #mod = lognorm(kx, 1.0, val, err)
        ky += mod
    ky /= np.float(len(values))
    
    return (kx, ky)


def wkde(values, errors, weights):

    minloc = np.argmin(values)
    minlim = np.clip(values[minloc] - 3*errors[minloc], 1e-6, np.inf)
    
    maxloc = np.argmax(values)
    maxlim = values[maxloc] + 3*errors[maxloc]
    
    #kx = np.linspace(minlim, maxlim, 1000)
    kx = np.logspace(np.log10(minlim), np.log10(maxlim), 1000)
    ky = np.zeros_like(kx)
        
    for i,val in enumerate(values):
        err = errors[i]
        w = weights[i]
        mod = w * gauss(kx, 1.0, val, err)            
        #mod = lognorm(kx, 1.0, val, err)
        ky += mod
    #ky /= np.float(len(values))
    
    return (kx, ky)


def wkde2D(xvalues, yvalues, xerrors, yerrors, weights, xlim=None, ylim=None, nstars=49545.):
    """Weighted KDE

    Calculate weighted KDE (wKDE) as described in the appendix of Fulton et al. (2017)

    Args:
        xvalues (array): array of x values (measurements)
        yvalues (array): array of y values (measurements)
        xerrors (array): uncertainties on x values
        yerrors (array): uncertainties on y values
        weights (array): array of weights for each measurement
        xlim (tuple): (optional) limits for x measurements, omit outside this range
        ylim (tuple): (optional) limits for y measurements, omit outside this range
        nstars (float): (optional) total number of stars observed in sample (default = 36075.0). Used to calculate
            absolute occurrence scale
    Returns:
        tuple: (x grid, y grid, z grid)

    """

    if xlim is None:
        xminloc = np.argmin(xvalues)
        xminlim = np.clip(xvalues[xminloc] - 10*xerrors[xminloc], 1e-6, np.inf)
        xmaxloc = np.argmax(xvalues)
        xmaxlim = xvalues[xmaxloc] + 10*xerrors[xmaxloc]
    else:
        xminlim, xmaxlim = xlim

    if ylim is None:
        yminloc = np.argmin(yvalues)
        yminlim = np.clip(yvalues[yminloc] - 10*yerrors[yminloc], 1e-6, np.inf)
        ymaxloc = np.argmax(yvalues)
        ymaxlim = yvalues[ymaxloc] + 10*yerrors[ymaxloc]
    else:
        yminlim, ymaxlim = ylim

    xspace = (np.log10(xmaxlim) - np.log10(xminlim)) / 512
    yspace = (np.log10(ymaxlim) - np.log10(yminlim)) / 512
    xi, yi = np.mgrid[np.log10(xminlim):np.log10(xmaxlim):xspace,
                      np.log10(yminlim):np.log10(ymaxlim):yspace]
    pos = np.dstack((xi, yi))
    
    z = np.zeros((len(xi), len(yi)))
    for i, val in enumerate(xvalues):
        x = xvalues[i]
        y = yvalues[i]
        # xe = np.clip((0.434*xerrors[i]/x)**2, 2.5e-5, np.inf)
        # ye = np.clip((0.434*yerrors[i]/y)**2, 2.5e-5, np.inf)
        xe = (0.434*xerrors[i]/x)**2
        ye = (0.434*yerrors[i]/y)**2

        w = weights[i]

        #print x, xe, y, ye, w
        mvn = scipy.stats.multivariate_normal([np.log10(x), np.log10(y)], [[xe, 0.0], [0.0, ye]])
        mod = mvn.pdf(pos)
        mod = w * (mod / np.max(mod))
        
        z += mod
    z /= np.float(nstars)
        
    return (xi, yi, z)
    
def draw_planets(binl, binr, values, num_draws=1000):
    N = np.array(np.round(values * num_draws), dtype=int)
    planets = []
    i = 0
    for l, r in zip(binl, binr):
        planets.append(np.random.uniform(l, r, size=N[i]))
        i += 1
        
    planets = np.hstack(planets)
        
    return planets
