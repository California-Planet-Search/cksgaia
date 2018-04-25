# Put plotting code here

from matplotlib import ticker
import matplotlib.lines as mlines

import os

import pylab as pl
import matplotlib
import numpy as np

import cksgaia.io
import cksgaia.plot.sample
from cksgaia.plot.config import modpath
from cksgaia.completeness import num_stars
from cksgaia.plot.config import *


def get_mass_samples():
    physmerge = cksgaia.io.load_table(cksgaia.plot.config.filtered_sample).query('giso_prad > 1.75 & giso_prad < 4')

    highcut = np.percentile(physmerge['giso_smass'], 67)
    lowcut = np.percentile(physmerge['giso_smass'], 33)

    high = physmerge.query('giso_smass > @highcut')
    medium = physmerge.query('giso_smass <= @highcut & giso_smass >= @lowcut')
    low = physmerge.query('giso_smass < @lowcut')

    annotations = ['$M_{\star} > %3.2f M_{\odot}$' % highcut,
                   '$%3.2f M_{\odot} \leq M_{\star} \leq %3.2f M_{\odot}$' % (highcut, lowcut),
                   '$M_{\star} < %3.2f M_{\odot}$' % lowcut]

    return (highcut, lowcut, high, medium, low, annotations)

def histfitplot(physmerge, bin_centers, N, e, mask, result, result2, completeness=False, plotmod=True,
                eloc=[10.0, 0.10], alpha=1.0, unc=True):
    xmod = np.logspace(-2, 1, 1000)

    masklim = (np.min(bin_centers[mask]), np.max(bin_centers[mask]))

    fig = pl.figure(figsize=(12, 8))

    pl.subplots_adjust(left=0.15)

    pl.step(bin_centers, N, color='k', where='mid', lw=4, alpha=alpha)
    # pl.bar(bin_centers, N, color='k', lw=4, alpha=alpha)
    pl.step(bin_centers[0:8], N[0:8], color='0.9', where='mid', lw=5, alpha=alpha)
    # pl.step(bin_centers[~mask][0:9], N[~mask][0:9], color='0.9', where='mid', lw=5, alpha=alpha)
    # pl.step(bin_centers[~mask][8:], N[~mask][8:], color='0.9', where='mid', lw=4)
    # pl.step(bin_centers[~mask][7:], NF[~mask][7:], color='0.9', where='mid', lw=3)

    e[e >= 1] = np.inf

    xpos = eloc[0]
    ypos = eloc[1]
    err1, err2 = cksgaia.misc.frac_err(physmerge, xpos, 'giso_prad')
    print "err+, err- = ", err1, err2

    _, caps, _ = pl.errorbar([xpos], [ypos], fmt='k.', xerr=[[err1], [err2]],
                             lw=2, capsize=6, mew=0, ms=0.1)
    for cap in caps:
        cap.set_markeredgewidth(2)
    pl.annotate("typical\nuncert.", xy=(xpos, ypos), xytext=(0, -50),
                xycoords="data", textcoords="offset points",
                horizontalalignment='center', fontsize=18)

    if unc:
        _, caps, _ = pl.errorbar(
            bin_centers[mask],
            N[mask],
            yerr=e[mask],
            marker='o',
            mew=0,
            drawstyle='steps-mid',
            linestyle='none',
            color='k',
            lw=2,
            capsize=6,
            alpha=alpha
        )

    for cap in caps:
        cap.set_markeredgewidth(2)

    # fit0 = cksrad.fitting.lognorm(xmod, result.params['amp'].value,
    #                          result.params['mu'].value,
    #                          result.params['sig'].value)

    fit1 = cksgaia.fitting.lognorm(xmod, result2.params['amp1'].value,
                                  result2.params['mu1'].value,
                                  result2.params['sig1'].value)
    fit2 = cksgaia.fitting.lognorm(xmod, result2.params['amp2'].value,
                                  result2.params['mu2'].value,
                                  result2.params['sig2'].value)

    # fit1_deconv = cksrad.fitting.lognorm(xmod, result2.params['amp1'].value,
    #                          result2.params['mu1'].value,
    #                          result2.params['sig1'].value, unc_limit=0.0)
    # fit1_deconv *= np.trapz(fit1)/np.trapz(fit1_deconv)
    # fit2_deconv = cksrad.fitting.lognorm(xmod, result2.params['amp2'].value,
    #                          result2.params['mu2'].value,
    #                          result2.params['sig2'].value, unc_limit=0.0)
    # fit2_deconv *= np.trapz(fit2)/np.trapz(fit2_deconv)

    # pl.plot(xmod, fit0, 'b--', lw=2, alpha=0.8)

    fit = cksgaia.fitting.splinefunc(xmod, result.params['n1'].value, result.params['n2'].value,
                                    result.params['n3'].value, result.params['n4'].value,
                                    result.params['n5'].value, result.params['n6'].value,
                                    result.params['n7'].value)
    fit_deconv = cksgaia.fitting.splinefunc(xmod, result.params['n1'].value, result.params['n2'].value,
                                           result.params['n3'].value, result.params['n4'].value,
                                           result.params['n5'].value, result.params['n6'].value,
                                           result.params['n7'].value, unc_limit=0)

    node_heights = [result.params['n%d' % (n + 1)].value for n in range(len(cksgaia.fitting.nodes))]

    if plotmod:
        c2 = (0, 146 / 255., 146 / 255.)
        c1 = (146 / 255., 0, 0)
        # pl.plot(xmod, fit1, '-', color=c1, lw=2, alpha=0.8)
        # pl.fill(xmod, fit1, color=c1, alpha=0.5)
        # pl.plot(xmod, fit2, '-', color=c2, lw=2, alpha=0.8)
        # pl.fill(xmod, fit2, color=c2, alpha=0.5)
        # pl.plot(xmod, fit1+fit2, 'r-', lw=3, alpha=0.8)
        # pl.plot(xmod, fit1_deconv+fit2_deconv, 'b-', lw=3, alpha=0.8)
        pl.plot(xmod, fit, '-', color=c1, lw=3)
        pl.plot(xmod, fit_deconv, '--', color=c2, lw=3)
        left = xmod <= masklim[0]
        right = xmod >= masklim[1]
        pl.plot(xmod[left], fit[left], '-', color='0.9', lw=3)
        pl.plot(xmod[right], fit_deconv[right], '--', color='0.9', lw=3)
        print cksgaia.fitting.nodes, node_heights
        pl.plot(10 ** cksgaia.fitting.nodes, node_heights, 'o', color='none', mew=2, ms=14, mec=c2)
    pl.semilogx()

    # pl.annotate("super-\nEarths", xy=(0.17, 0.07), xycoords='axes fraction')
    # pl.annotate("sub-\nNeptunes", xy=(0.45, 0.07), xycoords='axes fraction')

    pl.xlim(0.7, 20)
    pl.ylim(0, 0.12)
    if completeness:
        pl.ylabel('Number of Planets per Star\n(Orbital period < 100 days)')
    else:
        pl.ylabel('Number of Planets')
    pl.xlabel('Planet Size [Earth radii]')

    pl.xticks([.7, 1.0, 1.3, 1.8, 2.4, 3.5, 4.5, 6, 8.0, 12.0, 20.0])
    # pl.xticks(np.logspace(np.log10(0.7), np.log10(8), 6))

    ax = pl.gca()
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%0.1f'))

    return ax


def insol_hist():
    physmerge = cksgaia.io.load_table(cksgaia.plot.config.filtered_sample)

    # cx, cy = np.loadtxt('/Users/bfulton/code/cksrad/data/detectability_p1.txt', unpack=True)
    cx, cy = np.loadtxt(os.path.join(modpath, 'data/sensitivity_p25.txt'), unpack=True)
    # cx, cy = np.loadtxt('/Users/bfulton/code/cksrad/data/sensitivity_p50.txt', unpack=True)
    a = (physmerge['giso_smass'].max() * (cx / 365.) ** 2) ** (1 / 3.)
    sx = (physmerge['cks_steff'].max() / 5778) ** 4.0 * (physmerge['giso_srad'].max() / a) ** 2.0
    sx = np.append(sx, 10)
    cy = np.append(cy, 6)

    figure = pl.figure(figsize=(10, 12))
    nrow = 3
    ncol = 1
    plti = 1

    xticks = [.7, 1.0, 1.3, 1.75, 2.4, 3.5, 4.5, 6]

    pl.subplot(nrow, ncol, plti)
    pl.subplots_adjust(hspace=0, top=0.98, bottom=0.10, left=0.19)

    hot = physmerge.query('giso_insol > 200')
    medium = physmerge.query('giso_insol <= 200 & giso_insol >= 80')
    cool = physmerge.query('giso_insol < 50 & giso_insol > 10')

    lim = cy[np.argmin(np.abs(sx - 200))]
    print lim
    v = cksgaia.plot.sample.simplehist(hot, fill_valley=False, nbins=36, color='k', unc=True, aloc=(0.9, 0.8),
                                   annotate='$S_{\\rm inc} > 200 S_{\oplus}$', stacked=False, va_anno=False,
                                   weighted=True,
                                   nstars=num_stars, eloc=(4.5, 0.04), clim=lim)

    ax = pl.gca()
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)
    for label in ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    pl.xticks([])

    pl.ylabel('')
    pl.ylim(0, 0.06)
    pl.xlim(0.7, 8)
    pl.xticks([])
    plti += 1

    lim = cy[np.argmin(np.abs(sx - 80))]
    print lim
    pl.subplot(nrow, ncol, plti)
    v = cksgaia.plot.sample.simplehist(medium, fill_valley=False, nbins=36, color='k', unc=True, aloc=(0.9, 0.8),
                                   annotate='$80 S_{\oplus} \leq S_{\\rm inc} \leq 200 S_{\oplus}$', stacked=False,
                                   va_anno=False, weighted=True, nstars=num_stars, eloc=(4.5, 0.04), clim=lim)

    ax = pl.gca()
    # yticks = ax.yaxis.get_major_ticks()
    # yticks[0].label1.set_visible(False)
    for label in ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    ax.yaxis.get_ticklabels()[0].set_visible(False)
    ax.yaxis.get_ticklabels()[-1].set_visible(False)
    pl.xticks([])
    pl.ylabel('')
    pl.ylim(0, 0.06)
    pl.xlim(0.7, 8)
    pl.xticks([])
    plti += 1

    pl.ylabel('Number of Planets per Star')

    lim = cy[-2]
    print lim
    pl.subplot(nrow, ncol, plti)
    v = cksgaia.plot.sample.simplehist(cool, fill_valley=False, nbins=36, color='k', unc=True, aloc=(0.9, 0.85),
                                   annotate='$10 S_{\oplus} \leq S_{\\rm inc} \leq 50 S_{\oplus}$', stacked=False,
                                   va_anno=False, weighted=True, nstars=num_stars, eloc=(4.5, 0.04), clim=lim)

    ax = pl.gca()
    # yticks = ax.yaxis.get_major_ticks()
    # yticks[0].label1.set_visible(False)
    for label in ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    pl.xticks([])

    pl.ylabel('')
    pl.ylim(0, 0.06)
    pl.xlim(0.7, 8)
    pl.xticks(xticks)


def money_plot_fit():
    physmerge = cksgaia.io.load_table(cksgaia.plot.config.filtered_sample)

    rmask, rbin_centers, rN, re, result1, result2 = cksgaia.fitting.histfit(physmerge,
                                                                           completeness=False,
                                                                           verbose=False)

    mask, bin_centers, N, e, result1, result2 = cksgaia.fitting.histfit(physmerge, completeness=True, boot_errors=False)
    histfitplot(physmerge, bin_centers, N, e, mask, result1, result2, completeness=True)
    pl.grid(False)

    c2 = (0, 146 / 255., 146 / 255.)
    c1 = (146 / 255., 0, 0)

    # pl.axvline(1.0, color=c1, linestyle='dashed', lw=2)
    # pl.axvline(1.745, color=c1, linestyle='dashed', lw=2)
    pl.axvspan(1.0, 1.75, color=c1, alpha=0.1)

    # pl.axvline(1.760, color=c2, linestyle='dashed', lw=2)
    # pl.axvline(3.5, color=c2, linestyle='dashed', lw=2)
    pl.axvspan(1.75, 3.5, color=c2, alpha=0.1)

    pl.ylim(0, 0.125)


def money_plot_plain():
    physmerge = cksgaia.io.load_table(cksgaia.plot.config.filtered_sample)

    print len(physmerge), (physmerge.giso_srad_err1 / physmerge.giso_srad).median()

    rmask, rbin_centers, rN, re, result1, result2 = cksgaia.fitting.histfit(physmerge,
                                                                           completeness=False,
                                                                           verbose=False)

    mask, bin_centers, N, e, result1, result2 = cksgaia.fitting.histfit(physmerge, completeness=True, boot_errors=False)

    histfitplot(physmerge, bin_centers, N, e, mask, result1, result2, completeness=True, plotmod=False)
    pl.grid(False)

    rN = rN / (float(num_stars) * np.mean(physmerge['tr_prob']))
    pl.step(rbin_centers, rN, color='0.5', where='mid', lw=3, linestyle='dotted')

    c2 = (0, 146 / 255., 146 / 255.)
    c1 = (146 / 255., 0, 0)

    # pl.axvline(1.0, color=c1, linestyle='dashed', lw=2)
    # pl.axvline(1.745, color=c1, linestyle='dashed', lw=2)
    # pl.axvspan(1.0, 1.75, color=c1, alpha=0.1)

    # pl.axvline(1.760, color=c2, linestyle='dashed', lw=2)
    # pl.axvline(3.5, color=c2, linestyle='dashed', lw=2)
    # pl.axvspan(1.75, 3.5, color=c2, alpha=0.1)

    pl.ylim(0, 0.125)


def mass_cuts():
    physmerge = cksgaia.io.load_table(cksgaia.plot.config.filtered_sample)

    # cx, cy = np.loadtxt('/Users/bfulton/code/cksrad/data/detectability_p1.txt', unpack=True)
    cx, cy = np.loadtxt(os.path.join(modpath, 'data/sensitivity_p25.txt'), unpack=True)
    # cx, cy = np.loadtxt('/Users/bfulton/code/cksrad/data/sensitivity_p50.txt', unpack=True)
    a = (physmerge['giso_smass'].max() * (cx / 365.) ** 2) ** (1 / 3.)
    sx = (physmerge['cks_steff'].max() / 5778) ** 4.0 * (physmerge['giso_srad'].max() / a) ** 2.0
    sx = np.append(sx, 10)
    cy = np.append(cy, 6)

    figure = pl.figure(figsize=(10, 12))
    nrow = 3
    ncol = 1
    plti = 1

    xticks = [.7, 1.0, 1.3, 1.75, 2.4, 3.5, 4.5, 6]

    pl.subplot(nrow, ncol, plti)
    pl.subplots_adjust(hspace=0, top=0.98, bottom=0.10, left=0.19)

    highcut = np.percentile(physmerge['giso_smass'], 67)
    lowcut = np.percentile(physmerge['giso_smass'], 33)

    high = physmerge.query('giso_smass > @highcut')
    medium = physmerge.query('giso_smass <= @highcut & giso_smass >= @lowcut')
    low = physmerge.query('giso_smass < @lowcut')

    # lim = cy[np.argmin(np.abs(sx - 200))]
    lim = 1.14
    print lim
    v = cksgaia.plot.sample.simplehist(high, fill_valley=False, nbins=36, color='k', unc=True, aloc=(0.9, 0.8),
                                   annotate='$M_{\\star} > %3.2f M_{\odot}$' %highcut, stacked=False, va_anno=False,
                                   weighted=True,
                                   nstars=num_stars, eloc=(4.5, 0.04), clim=lim)

    ax = pl.gca()
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)
    for label in ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    pl.xticks([])

    pl.ylabel('')
    pl.ylim(0, 0.06)
    pl.xlim(0.7, 8)
    pl.xticks([])
    plti += 1

    # lim = cy[np.argmin(np.abs(sx - 80))]
    lim = 1.14
    print lim
    pl.subplot(nrow, ncol, plti)
    v = cksgaia.plot.sample.simplehist(medium, fill_valley=False, nbins=36, color='k', unc=True, aloc=(0.9, 0.8),
                                   annotate='$%3.2f M_{\odot} \leq M_{\\star} \leq %3.2f M_{\odot}$' % (highcut, lowcut),
                                   stacked=False, va_anno=False, weighted=True, nstars=num_stars,
                                   eloc=(4.5, 0.04), clim=lim)

    ax = pl.gca()
    # yticks = ax.yaxis.get_major_ticks()
    # yticks[0].label1.set_visible(False)
    for label in ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    ax.yaxis.get_ticklabels()[0].set_visible(False)
    ax.yaxis.get_ticklabels()[-1].set_visible(False)
    pl.xticks([])
    pl.ylabel('')
    pl.ylim(0, 0.06)
    pl.xlim(0.7, 8)
    pl.xticks([])
    plti += 1

    pl.ylabel('Number of Planets per Star')

    lim = cy[-2]
    lim = 1.14
    # print lim
    pl.subplot(nrow, ncol, plti)
    v = cksgaia.plot.sample.simplehist(low, fill_valley=False, nbins=36, color='k', unc=True, aloc=(0.9, 0.85),
                                   annotate='$M_{\star} < %3.2f M_{\odot}$' % lowcut, stacked=False,
                                   va_anno=False, weighted=True, nstars=num_stars, eloc=(4.5, 0.04), clim=lim)

    ax = pl.gca()
    # yticks = ax.yaxis.get_major_ticks()
    # yticks[0].label1.set_visible(False)
    for label in ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    pl.xticks([])

    pl.ylabel('')
    pl.ylim(0, 0.06)
    pl.xlim(0.7, 8)
    pl.xticks(xticks)

def desert_edge():
    physmerge = cksgaia.io.load_table(cksgaia.plot.config.filtered_sample).query('giso_prad > 1.75 & giso_prad < 4')

    aloc = (0.1, 0.85)

    # cx, cy = np.loadtxt('/Users/bfulton/code/cksrad/data/detectability_p1.txt', unpack=True)
    cx, cy = np.loadtxt(os.path.join(modpath, 'data/sensitivity_p25.txt'), unpack=True)
    # cx, cy = np.loadtxt('/Users/bfulton/code/cksrad/data/sensitivity_p50.txt', unpack=True)
    a = (physmerge['giso_smass'].max() * (cx / 365.) ** 2) ** (1 / 3.)
    sx = (physmerge['cks_steff'].max() / 5778) ** 4.0 * (physmerge['giso_srad'].max() / a) ** 2.0
    sx = np.append(sx, 10)
    cy = np.append(cy, 6)

    figure = pl.figure(figsize=(10.1, 12))
    nrow = 3
    ncol = 1
    plti = 1

    xticks = [10000, 3000, 1000, 300, 100, 30, 10]

    pl.subplot(nrow, ncol, plti)
    pl.subplots_adjust(hspace=0, top=0.98, bottom=0.10, left=0.19)

    highcut = np.percentile(physmerge['giso_smass'], 67)
    lowcut = np.percentile(physmerge['giso_smass'], 33)

    high = physmerge.query('giso_smass > @highcut')
    medium = physmerge.query('giso_smass <= @highcut & giso_smass >= @lowcut')
    low = physmerge.query('giso_smass < @lowcut')

    # lim = cy[np.argmin(np.abs(sx - 200))]
    lim = 1.14
    print lim

    insolbins = np.logspace(np.log10(10), np.log10(10000), 20)

    cut = high
    N, edges = np.histogram(cut['giso_insol'].dropna().values, bins=insolbins, weights=cut['weight'])
    Nd, edges = np.histogram(cut['giso_insol'].dropna().values, bins=insolbins)
    centers = 0.5 * (edges[1:] + edges[:-1])
    pl.step(centers, N / num_stars, lw=3, color='k', where='mid')
    pl.annotate('$M_{\\star} > %3.2f M_{\odot}$' %highcut, xy=aloc, xycoords='axes fraction', fontsize=20,
                horizontalalignment='left')

    err = (np.sqrt(Nd) * (N / Nd)) / num_stars
    err[np.isnan(err)] = 0
    _, caps, _ = pl.errorbar(centers, N / num_stars, fmt='.', yerr=err, lw=2, capsize=6, color='k')
    for cap in caps:
        cap.set_markeredgewidth(2)

    pl.semilogx()


    ax = pl.gca()
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%0.0f'))
    ax.xaxis.set_ticks(np.logspace(1, 4, 8))

    yticks = ax.yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)
    for label in ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    for label in ax.xaxis.get_ticklabels():
        label.set_visible(False)
    pl.xticks(xticks)


    pl.ylabel('')
    #pl.ylim(0, 0.06)
    pl.xlim(10000, 10)
    #pl.xticks([])
    plti += 1

    # lim = cy[np.argmin(np.abs(sx - 80))]
    lim = 1.14
    print lim
    pl.subplot(nrow, ncol, plti)

    cut = medium
    N, edges = np.histogram(cut['giso_insol'].dropna().values, bins=insolbins, weights=cut['weight'])
    Nd, edges = np.histogram(cut['giso_insol'].dropna().values, bins=insolbins)
    centers = 0.5 * (edges[1:] + edges[:-1])
    pl.step(centers, N / num_stars, lw=3, color='k', where='mid')
    pl.annotate('$%3.2f M_{\odot} \leq M_{\\star} \leq %3.2f M_{\odot}$' % (highcut, lowcut), xy=aloc, xycoords='axes fraction', fontsize=20,
                horizontalalignment='left')

    err = (np.sqrt(Nd) * (N / Nd)) / num_stars
    err[np.isnan(err)] = 0
    _, caps, _ = pl.errorbar(centers, N / num_stars, fmt='.', yerr=err, lw=2, capsize=6, color='k')
    for cap in caps:
        cap.set_markeredgewidth(2)

    pl.semilogx()

    ax = pl.gca()
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%0.0f'))
    ax.xaxis.set_ticks(np.logspace(1, 4, 8))

    # yticks = ax.yaxis.get_major_ticks()
    # yticks[0].label1.set_visible(False)
    for label in ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    for label in ax.xaxis.get_ticklabels():
        label.set_visible(False)
    ax.yaxis.get_ticklabels()[0].set_visible(False)
    ax.yaxis.get_ticklabels()[-1].set_visible(False)
    pl.xticks(xticks)
    pl.ylabel('')
    #pl.ylim(0, 0.06)
    pl.xlim(10000, 10)
    pl.xticks(xticks)
    plti += 1

    pl.ylabel('Number of Planets per Star')

    lim = cy[-2]
    lim = 1.14
    # print lim
    pl.subplot(nrow, ncol, plti)

    cut = low
    N, edges = np.histogram(cut['giso_insol'].dropna().values, bins=insolbins, weights=cut['weight'])
    Nd, edges = np.histogram(cut['giso_insol'].dropna().values, bins=insolbins)
    centers = 0.5 * (edges[1:] + edges[:-1])
    pl.step(centers, N / num_stars, lw=3, color='k', where='mid')
    pl.annotate('$M_{\star} < %3.2f M_{\odot}$' % lowcut, xy=aloc, xycoords='axes fraction', fontsize=20,
                horizontalalignment='left')

    err = (np.sqrt(Nd) * (N / Nd)) / num_stars
    err[np.isnan(err)] = 0
    _, caps, _ = pl.errorbar(centers, N / num_stars, fmt='.', yerr=err, lw=2, capsize=6, color='k')
    for cap in caps:
        cap.set_markeredgewidth(2)

    pl.semilogx()

    ax = pl.gca()
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%0.0f'))
    ax.xaxis.set_ticks(np.logspace(1, 4, 8))

    # yticks = ax.yaxis.get_major_ticks()
    # yticks[0].label1.set_visible(False)
    for label in ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    pl.xticks([])

    pl.ylabel('')
    #pl.ylim(0, 0.06)
    pl.xlim(10000, 10)
    pl.xticks(xticks)

    pl.xlabel('Stellar light intensity relative to Earth')


def desert_edge_cum():

    def _cumdist(sample, annotation='', color='k'):
        order = np.argsort(sample['giso_insol'].values)[::-1]
        w = np.cumsum(sample['weight'].values[order]) / num_stars
        n = np.array(range(len(sample['weight'].values))) / float(len(sample))
        f = sample['giso_insol'].values[order]

        pl.step(f, n, 'k-', lw=3, linestyle='dashed', alpha=0.25, color=color, label=None)
        pl.step(f, w / np.max(w), 'k-', lw=3, color=color)

        aloc = (0.1, 0.85)
        #pl.annotate(annotation, xy=aloc, xycoords='axes fraction', fontsize=20,
        #            horizontalalignment='left')

        pl.semilogx()

        pl.xlim(3000, 10)
        pl.ylim(0, 1)

        ax = pl.gca()
        ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%0.0f'))
        ax.xaxis.set_ticks([10, 30, 100, 300, 1000, 3000])

        pl.xlabel('Stellar light intensity relative to Earth')
        pl.ylabel('cumulative fraction of planet detections')

    highcut, lowcut, high, medium, low, annotations = get_mass_samples()

    colors = ['blue', 'green', 'red']

    fig = pl.figure(figsize=(12, 8))
    handles = []
    for i, sample in enumerate([high, medium, low]):
        sample = sample.query('giso_prad > 1.75 & giso_prad < 4')

        _cumdist(sample, color=colors[i])
        handles.append(mlines.Line2D([], [], color=colors[i], lw=3,
                                     label=annotations[i]))

    pl.legend(handles=handles, fontsize=14, loc='best')

    pl.annotate('$1.75 < R_p < 4 R_{\oplus}$', xy=(0.03, 0.75),
                xycoords='axes fraction', fontsize=16)


