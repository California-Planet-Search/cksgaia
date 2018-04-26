
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt

import cksgaia.io
import cksgaia.calc
import cksgaia.completeness

from cksgaia.config import *
from cksgaia.plot.config import *

matplotlib.rcParams['figure.figsize'] = (13, 8)

def contour_plot_kde(physmerge, xcol, ycol, xlim, ylim, ylog=True, pltxlim=None, pltylim=None, epos=[3000, 5],
                     cont=True, nodata=False, weighted=False, nstars=36075., eaoff=(0, -70), clabel=None,
                     vlims=(0.0, 0.05), kwidth=None):
    """Plot contour plots

    Make the contour plots associated with Figures 8-10 in Fulton et al. (2017)

    Args:
        physmerge (DataFrame): Pandas DataFrame with the input planet catalogue
        xcol (string): name of x-axis column to plot
        ycol (string): name of y-axis column to plot
        xlim (tuple): x axis limits (min, max) for the wKDE calculation, should be slightly outside the pltxlim parameter
        ylim (tuple): y axis limits (min, max) for the wKDE calculation, should be slightly outside the pltylim parameter
        ylog (bool): (optional) log scale for y-axis?
        pltxlim (tuple): (optional) x-axis limit to display on plot (if different from xlim)
        pltylim (tuple): (optional) y-axis limit to display on plot (if different from ylim)
        epos (tuple): (optional) position to plot the typical uncertainty error bar
        cont (bool): (optional) plot the contours (default = True)
        nodata (bool): (optional) omit the datapoints from the plot (default = False)
        weighted (bool): (optional) use the 'weight' column to calculate a weighted KDE (default = False)
        nstars (float): (optional) total number of stars observed in sample (default = 36075.0). Used to calculate
            absolute occurrence scale
        eaoff (tuple): (optional) text offset for label of the typical uncertainty error bar (default = (0, -70))
        clabel (string): (optional) string to label color bar scale
        vlims (tuple): (optional) colorscale limits (default = (0.0, 0.05))
        kwidth (tuple): (optional) kernel width in units of (x,y)

    Returns:
        tuple: (axes object, 10**(grid x values), 10**(grid y values), grid z values)

    """

    # fig = pl.figure(figsize=(13, 8))

    crop = physmerge[(physmerge[ycol] < ylim[1]) & (physmerge[ycol] > ylim[0]) &
                     (physmerge[xcol] < xlim[1]) & (physmerge[xcol] < xlim[1])]

    nplanets = float(len(crop))

    if weighted:
        weights = crop['weight'].values
    else:
        weights = np.ones_like(crop[xcol].values)

    xerr1 = crop[xcol + '_err1'].values
    yerr1 = crop[ycol + '_err1'].values
    if kwidth != None:
        print kwidth
        xerr1 = np.zeros_like(xerr1) + kwidth[0] * crop[xcol].values
        yerr1 = np.zeros_like(yerr1) + kwidth[1] * crop[ycol].values

    xi, yi, zi = cksgaia.fitting.wkde2D(crop[xcol].values, crop[ycol].values,
                                       xerr1, yerr1,
                                       weights, xlim=xlim, ylim=ylim, nstars=nstars)

    # cmap = plt.cm.inferno_r
    # cmap = plt.cm.gray_r
    # cmap = plt.cm.bone_r
    # cmap = plt.cm.gist_heat_r
    # cmap = plt.cm.magma_r
    cmap = plt.cm.afmhot_r
    # cmap = plt.cm.Blues_r

    vmin, vmax = vlims
    levels = np.arange(vmin, vmax + 0.0001, 0.005)
    if len(levels) <= 7:
        levels = np.linspace(vmin, vmax + 1e-4, 11)

    if cont:
        CS = plt.contourf(10 ** xi, 10 ** yi, zi, levels=levels, cmap=cmap, vmin=vmin, vmax=vmax)
        for c in CS.collections:
            c.set_edgecolor("face")

    # pl.plot(big_sdist, big_pdist, 'k.', color='0.3', ms=1.0)
    if not nodata:
        pl.plot(crop[xcol], crop[ycol], 'ko', ms=6, markeredgewidth=1, markeredgecolor='w')

    xe, ye = epos
    xerr1, xerr2 = cksgaia.misc.frac_err(physmerge, xe, xcol)
    yerr1, yerr2 = cksgaia.misc.frac_err(physmerge, ye, ycol)

    _, caps, _ = pl.errorbar(xe, ye, yerr=[[yerr1], [yerr2]], xerr=[[xerr1], [xerr2]],
                             markeredgewidth=0, lw=2, capsize=6, markeredgecolor='w', fmt='k.', ms=0.1)
    for cap in caps:
        cap.set_markeredgewidth(2)
    pl.annotate("typical\nuncert.", xy=(xe, ye), xytext=eaoff,
                xycoords="data", textcoords="offset points",
                horizontalalignment='center', fontsize=18)

    # pl.semilogx()
    if ylog:
        pl.loglog()
        pl.xlim(pltxlim)
        pl.ylim(pltylim)
    else:
        pl.semilogx()
        pl.xlim(pltxlim)
        pl.ylim(pltylim)

    yticks = np.array([0.3, 1, 2, 4, 6, 10, 30])
    xticks = np.array([0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 3e4])
    yt = yticks[(yticks >= ylim[0]) & (yticks <= ylim[1])]
    xt = xticks[(xticks >= xlim[0]) & (xticks <= xlim[1])]

    ax = pl.gca()
    if ylog:
        pl.yticks(yt)
        pl.xticks(xt)
    else:
        pl.xticks(xt)

    pl.xlim(pltxlim)
    pl.ylim(pltylim)

    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.0f'))
    if ylog:
        ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.0f'))
    else:
        ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))

    pl.grid(False)

    if clabel is not None:
        clabel = clabel
    elif weighted:
        clabel = 'Relative Occurrence'
    else:
        clabel = 'Relative Density of Planets'
    cmap = pl.colorbar(pad=0, label=clabel)

    return (ax, 10 ** xi, 10 ** yi, zi)


def period_contour_q16():
    physmerge = cksgaia.io.load_table(cksgaia.plot.config.filtered_sample)

    ax, xi, yi, zi = contour_plot_kde(physmerge, 'koi_period', 'koi_prad', xlim=[0.4, 1000.0],
                                                      ylim=[0.5, 20], ylog=True,
                                                      pltxlim=[0.7, 100.0], pltylim=[1.0, 10], epos=[1.2, 6.0],
                                                      eaoff=(0, -100), weighted=True,
                                                      vlims=(0.0, 0.07))

    # pl.xticks([0.6,0.8,1.0,1.2,1.5])
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.0f'))

    pl.xlabel('Orbital period [days]')
    pl.ylabel('Planet Size [Earth radii]')

    pl.title('Q16')


def period_contour_cks(sample=None, kwidth=(0.40, 0.05), vlims=(0.0, 0.05),
                       ylimits=(1.0, 10.0), clim=None):
    if sample is None:
        physmerge = cksgaia.io.load_table(cksgaia.plot.config.filtered_sample)
    else:
        physmerge = sample

    wper, wsens = np.genfromtxt(os.path.join(modpath, 'data/detectability_p1.txt'), unpack=True)

    ax, xi, yi, zi = contour_plot_kde(physmerge, 'koi_period', 'giso_prad', xlim=[0.4, 1000.0],
                                                      ylim=[0.5, 20], ylog=True,
                                                      pltxlim=[0.7, 100.0], pltylim=ylimits, epos=[1.2, 8.0],
                                                      weighted=True, kwidth=kwidth, vlims=vlims)

    # pl.plot(wper, wsens, 'k-', color='0.6', linestyle='dashed', lw=3)

    # pl.xticks([0.6,0.8,1.0,1.2,1.5])
    pl.yticks([1.0, 1.5, 2.4, 3.5])
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.0f'))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))

    if clim is not None:
        cx, cy = clim
    else:
        cx, cy = np.loadtxt(os.path.join(DATADIR, 'sensitivity_p25.txt'), unpack=True)

    pl.fill_between(cx, cy, y2=0.1, color='0.2', zorder=10, alpha=0.8, hatch='\\\\')
    # pl.annotate('      low\ncompleteness', xy=(30, 1.03),
    #             xycoords='data', color='0.2', fontsize=afs - 2)

    # for v in [1000, 300, 100, 30, 10]:
    #     xs = xi[:,0].flatten()
    #     print xs.shape
    #     k = cksgaia.calc.kdeslice(xs, v, zi)
    #
    #     print k.shape
    #     np.savetxt(os.path.join(DATADIR, 'kde_slices_{}.txt'.format(v)),
    #                np.transpose((xi, k)))

    pl.xlabel('Orbital period [days]')
    pl.ylabel('Planet Size [Earth radii]')
    # pl.title('CKS')


def insol_contour_anno():
    physmerge = cksgaia.io.load_table(cksgaia.plot.config.filtered_sample)

    ax, xi, yi, zi_iso = contour_plot_kde(physmerge, 'giso_insol', 'giso_prad', xlim=[3, 30000],
                                                          ylim=[0.5, 10], ylog=True,
                                                          pltxlim=[10, 3000], pltylim=[1, 4], epos=[2000, 1.3],
                                                          nodata=True, weighted=True)

    matplotlib.rcParams['font.size'] = 24

    ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
    pl.yticks([1.0, 1.5, 2.4, 3.5])

    # pl.plot([100.0,500], [4,1.7], color='b', linestyle='dashed', lw=4)
    pl.plot([650.0, 3000], [2.2, 2.2], color='b', linestyle='dashed', lw=4)
    pl.plot([650.0, 3000], [3.8, 3.8], color='b', linestyle='dashed', lw=4)
    pl.plot([650.0, 650], [2.2, 3.8], color='b', linestyle='dashed', lw=4)
    # pl.annotate('photoevaporation desert', xy=(600,3.5), rotation=59, xycoords='data', color='b', fontsize=afs+1)
    pl.annotate('photoevap.\ndesert', xy=(2500, 3.0), rotation=0, xycoords='data', color='b', fontsize=afs + 1)
    # pl.plot([30, 1000], [1.9,1.9], 'w--', lw=4)

    # cx, cy = np.loadtxt('/Users/bfulton/code/cksrad/data/detectability_p1.txt', unpack=True)
    cx, cy = np.loadtxt(os.path.join(modpath, 'data/sensitivity_p25.txt'), unpack=True)
    # cx, cy = np.loadtxt('/Users/bfulton/code/cksrad/data/sensitivity_p50.txt', unpack=True)
    a = (physmerge['giso_smass'].max() * (cx / 365.) ** 2) ** (1 / 3.)
    sx = (physmerge['cks_steff'].max() / 5778) ** 4.0 * (physmerge['giso_srad'].max() / a) ** 2.0
    sx = np.append(sx, 10)
    cy = np.append(cy, 6)

    pl.fill_between(sx, cy, y2=0.1, color='0.5', zorder=10, alpha=0.5, hatch='\\\\')
    pl.annotate('      low\ncompleteness', xy=(30, 1.03), xycoords='data', color='0.2', fontsize=afs - 2)
    # pl.fill_between([10.0,180], [3.0,1.0], y2=0.1, color='0.5', zorder=10, alpha=0.5, hatch='\\\\')
    # pl.annotate('low completeness', xy=(70, 1.1), xycoords='data', color='0.2', fontsize=afs)

    vf = np.logspace(np.log10(10), np.log10(325), 100)
    vr_evap = vf ** 0.11
    vfi = np.logspace(np.log10(10), np.log10(500), 100)
    vr_impacts = vfi ** (-0.065) + 0.9

    pl.plot(vf, vr_evap, '--', color='purple', lw=4)
    pl.annotate('atmospheric loss', xy=(280, 1.90), xycoords='data', rotation=-18, color='purple', fontsize=afs)
    pl.plot(vfi, vr_impacts, linestyle='dotted', color='k', lw=4)
    pl.annotate('gas-poor', xy=(30, 1.63), xycoords='data', rotation=3.2, color='k', fontsize=afs)

    pl.xlabel('Stellar light intensity relative to Earth')
    pl.ylabel('Planet Size [Earth radii]')

    pl.xlim(pl.xlim()[::-1])


def insol_contour_data(sample=None, vlims=None, kwidth=(0.4, 0.05), clims=None):
    if sample is None:
        physmerge = cksgaia.io.load_table(cksgaia.plot.config.filtered_sample)
    else:
        physmerge = sample

    # cx, cy = np.loadtxt('/Users/bfulton/code/cksrad/data/detectability_p1.txt', unpack=True)
    if clims is None:
        cx, cy = np.loadtxt(os.path.join(modpath, 'data/sensitivity_p25.txt'), unpack=True)
    else:
        cx, cy = clims
    # cx, cy = np.loadtxt('/Users/bfulton/code/cksrad/data/sensitivity_p50.txt', unpack=True)
    a = (physmerge['giso_smass'].max() * (cx / 365.) ** 2) ** (1 / 3.)
    sx = (physmerge['cks_steff'].max() / 5778) ** 4.0 * (physmerge['giso_srad'].max() / a) ** 2.0
    sx = np.append(sx, 10)
    cy = np.append(cy, 6)

    if vlims is None:
        vlims = (0.0, 0.035)
    else:
        vlims = vlims

    ax, xi, yi, zi_iso = contour_plot_kde(physmerge, 'giso_insol', 'giso_prad', xlim=[3, 30000],
                                                          ylim=[0.5, 20], ylog=True,
                                                          pltxlim=[10, 3000], pltylim=[1, 4], epos=[1800, 3.0],
                                                          weighted=True, vlims=vlims, kwidth=kwidth)

    pl.fill_between(sx, cy, y2=0.1, color='0.2', zorder=10, alpha=0.8, hatch='\\\\')
    # pl.annotate('      low\ncompleteness', xy=(50, 1.03), xycoords='data', color='0.2', fontsize=afs - 2)

    pl.xlabel('Stellar light intensity relative to Earth')
    pl.ylabel('Planet Size [Earth radii]')

    pl.xlim(pl.xlim()[::-1])

    ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
    pl.yticks([])
    pl.yticks([1.0, 1.5, 2.4, 3.5])


def srad_contour():
    physmerge = cksgaia.io.load_table(cksgaia.plot.config.filtered_sample)

    ax, xi, yi, zi = contour_plot_kde(physmerge, 'giso_srad', 'giso_prad', xlim=[0.4, 3.0],
                                                      ylim=[0.5, 20], ylog=True,
                                                      pltxlim=[0.6, 2.3], pltylim=[1.0, 5], epos=[0.7, 4.0],# nbins=30,
                                                      weighted=False,
                                                      clabel="Relative Density of Planets", vlims=(0.0, 0.0035))
    pl.xlabel('Stellar Radius [Solar radii]')
    pl.ylabel('Planet Size [Earth radii]')

    pl.xticks([0.6, 0.8, 1.0, 1.2, 1.5, 2.0])
    pl.yticks([1, 2, 3, 4])

    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))


def contour_masscuts():
    physmerge = cksgaia.io.load_table(cksgaia.plot.config.filtered_sample)
    kicselect = cksgaia.io.load_table('kic-filtered')
    kicselect = cksgaia.completeness.fit_cdpp(kicselect)

    highcut, lowcut, _, _, _, annotations = cksgaia.plot.occur.get_mass_samples()

    high = physmerge.query('giso_smass > @highcut')
    medium = physmerge.query('giso_smass <= @highcut & giso_smass >= @lowcut')
    low = physmerge.query('giso_smass < @lowcut')

    fig = pl.figure(1, figsize=(36, 8))

    pl.subplot(1, 3, 1)

    pl.subplots_adjust(left=0.05, right=0.95)

    vlimits = [(0.0, 0.03), (0.0, 0.02), (0.0, 0.03)]

    for i, sample in enumerate([high, medium, low]):

        pl.subplot(1, 3, i+1)

        if sample is high:
            kicsample = kicselect.query('kic_smass > @highcut')
        elif sample is medium:
            kicsample = kicselect.query('kic_smass <= @highcut & kic_smass >= @lowcut')
        elif sample is low:
            kicsample = kicselect.query('kic_smass < @lowcut')

        sample = cksgaia.completeness.get_weights(sample, kicsample)

        # get 25% completeness limits
        pl.figure(2)
        cx, cy = cksgaia.completeness.get_sensitivity_contour(kicsample, 0.25)

        pl.figure(1)
        insol_contour_data(sample=sample, vlims=vlimits[i], kwidth=(0.75, 0.05))
        pl.title(annotations[i], fontsize=afs + 6)
        pl.grid(lw=2, alpha=0.5)


def period_contour_masscuts():
    physmerge = cksgaia.io.load_table(cksgaia.plot.config.filtered_sample)
    kicselect = cksgaia.io.load_table('kic-filtered')
    kicselect = cksgaia.completeness.fit_cdpp(kicselect)

    highcut, lowcut, _, _, _, annotations = cksgaia.plot.occur.get_mass_samples()

    high = physmerge.query('giso_smass > @highcut')
    medium = physmerge.query('giso_smass <= @highcut & giso_smass >= @lowcut')
    low = physmerge.query('giso_smass < @lowcut')

    fig = pl.figure(1, figsize=(36, 8))

    pl.subplot(1, 3, 1)
    pl.subplots_adjust(left=0.05, right=0.95)

    fig2 = pl.figure(2)

    vlimits = [(0.0, 0.04), (0.0, 0.023), (0.0, 0.04)]

    for i, sample in enumerate([high, medium, low]):

        if sample is high:
            kicsample = kicselect.query('kic_smass > @highcut')
        elif sample is medium:
            kicsample = kicselect.query('kic_smass <= @highcut & kic_smass >= @lowcut')
        elif sample is low:
            kicsample = kicselect.query('kic_smass < @lowcut')

        sample = cksgaia.completeness.get_weights(sample, kicsample)

        # get 25% completeness limits
        pl.figure(2)
        cx, cy = cksgaia.completeness.get_sensitivity_contour(kicsample, 0.25)

        pl.figure(1)
        pl.subplot(1, 3, i+1)
        period_contour_cks(sample=sample, vlims=vlimits[i],
                           kwidth=(0.75, 0.05), ylimits=(1.0, 4.0), clim=(cx, cy))
        pl.title(annotations[i], fontsize=afs + 6)
        pl.grid(lw=2, alpha=0.5)
