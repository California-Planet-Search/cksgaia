#!/usr/bin/env python

import numpy as np
import pylab as pl
import matplotlib
import pp
import time
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.stats import lognorm

import simulations
import cksgaia.fitting

nsim = 100
npl = 65000

matplotlib.rcParams['font.size'] = 24
matplotlib.rcParams['figure.figsize'] = (12,8)



try:
    server.destroy()
except:
    pass

server = pp.Server(ncpus=100)

time.sleep(5)
print server.get_active_nodes(), np.sum(server.get_active_nodes().values())

jobs = []
for i in range(nsim):
    jobs.append(server.submit(simulations.simulated_survey, args=(npl,)))
    #out = cksrad.simulations.simulated_survey(1000)

true_ratios = []
final_ratios = []
whistograms = []
rhistograms = []
ndetections = []
det_histograms = []
dfs = []
outdicts = []
for j in jobs:
    out = j()

    true_ratios.append(out['true_ratio'])
    final_ratios.append(out['final_ratio'])
    whistograms.append(out['weighted_hist_norm'])
    rhistograms.append(out['observed_hist_norm'])
    ndetections.append(out['n_detections'])
    dfs.append(out['sim_data'])
    det_histograms.append(out['observed_hist'])
    outdicts.append(out)

server.destroy()


hkx = np.logspace(np.log10(0.5), np.log10(8), 1000)

# Make some fake planets
r1frac = 0.5
r2frac = 1-r1frac
nstars = 5e4

# radii limits to define two peaks
r1, r2, r3, r4 = 1.1, 1.7, 2.0, 3.0

# broad distributions
# simulated_periods = np.random.lognormal(3.5, 1.0, size=int(np.round(npl*r1frac)))
# simulated_periods = np.append(simulated_periods, np.random.lognormal(2.5, 0.9, size=int(np.round(npl*r2frac))))
# simulated_radii = np.random.lognormal(1.05, 0.32, size=int(np.round(npl*r1frac)))
# simulated_radii = np.append(simulated_radii, np.random.lognormal(0.0, 0.20, size=int(np.round(npl*r2frac))))
# ly1 = cksgaia.fitting.lognorm(hkx, 1.0, 2.858, 0.32)
# ly2 = cksgaia.fitting.lognorm(hkx, 1.3, 0.819, 0.28)


# narrow distributions
# simulated_periods = np.random.lognormal(3.5, 0.7, size=int(np.round(npl*r1frac)))
# simulated_periods = np.append(simulated_periods, np.random.lognormal(2.5, 0.9, size=int(np.round(npl*r2frac))))
# simulated_radii = np.random.lognormal(0.867, 0.17, size=int(np.round(npl*r1frac)))
# simulated_radii = np.append(simulated_radii, np.random.lognormal(0.26, 0.12, size=int(np.round(npl*r2frac))))
ly1 = cksgaia.fitting.lognorm(hkx, 1.0, 2.38, 0.15)
ly2 = cksgaia.fitting.lognorm(hkx, 1.41, 1.30, 0.09)
ly1 = cksgaia.fitting.gauss(np.log10(hkx), 1.0, np.log10(2.38), np.log10(1.185))
ly2 = cksgaia.fitting.gauss(np.log10(hkx), 1.41, np.log10(1.30), np.log10(1.13))


# uniform distributions
simulated_periods = 10**np.random.uniform(np.log(1), np.log(200), size=int(npl))
simulated_radii = 10**(np.random.uniform(np.log(0.5), np.log(8), size=int(npl)))



random_order = np.array(range(npl))
np.random.shuffle(random_order)
simulated_periods = simulated_periods[random_order]
simulated_radii = simulated_radii[random_order]

Redges = np.logspace(np.log10(0.5), np.log10(20), 36)
true_dist,_ = np.histogram(simulated_radii, bins=Redges)
true_dist = np.array(true_dist, dtype=float) / nstars

ly = ly1 + ly2

print np.trapz(ly1, np.log10(hkx)), np.trapz(ly2, np.log10(hkx))


rcen = out['radius_bin_centers']

hstack = np.vstack(whistograms).transpose()
rstack = np.vstack(det_histograms).transpose()
med_det_norm = np.vstack(rhistograms).transpose().mean(axis=1)
print hstack.shape

rerr = hstack.std(axis=1)
med_hist = hstack.mean(axis=1)
med_det = rstack.mean(axis=1)
print med_hist.shape, rerr.shape


#print np.trapz(med_hist, rcen), np.trapz(ly, hkx)
tly = ly / (np.trapz(ly, hkx) / np.trapz(med_hist, rcen))


pl.step(rcen, med_det_norm, lw=4, color='r', where='mid', linestyle='--',
        label='simulated detections')
pl.step(rcen, hstack, lw=1, color='0.5', where='mid', linestyle='-', alpha=0.2)
pl.step(rcen, med_hist, lw=4, color='k', where='mid', linestyle='-',
            label='recovered distribution')
# pl.step(rcen, true_dist, lw=3, color='b', where='mid', linestyle='-',
#              label='true distribution')
pl.step(hkx, tly, lw=3, color='b', where='mid', linestyle='-',
            label='true distribution')

(_, caps, _) = pl.errorbar(rcen, med_hist, yerr=rerr, fmt='k.', lw=3, capsize=6)
for cap in caps:
    cap.set_markeredgewidth(3)
#pl.errorbar(Rcen, whistn, yerr=werr, fmt='k.')

sfac = 1/out['sim_data'].tr_prob.mean()
poisson_errors = np.sqrt(med_det) / 5e4 * sfac
poisson_errors *= med_hist/med_det_norm
err_ratio = poisson_errors/rerr
print repr(1/err_ratio), 1/np.nanmean(err_ratio)
pl.errorbar(rcen, med_hist, yerr=poisson_errors, fmt='r.', lw=2, capsize=6)


# pl.axvline(1.1)
# pl.axvline(1.7)
# pl.axvline(2.0)
# pl.axvline(3.0)
# print np.log10(3.) - np.log10(2)
# print np.log10(1.7) - np.log10(1.1)

pl.xlim(0.5,6)
pl.ylim(0,0.17)
pl.semilogx()
ax = pl.gca()
ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))

pl.xticks(np.array([0.5,0.7,1.0,1.5,2.2,3.2,5,8.0]))

pl.ylabel('Number of Planets per Star')
pl.xlabel('Planet Radius [Earth radii]')

pl.legend(loc='upper right')

print np.sum(med_hist)

pl.grid(False)

pl.savefig('/home/bfulton/Dropbox/plots/simulated_surveys_hist.pdf')


