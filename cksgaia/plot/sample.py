import pylab as pl
import numpy as np
import matplotlib

import cksgaia.io
import cksgaia.misc
import cksgaia.errors
from cksgaia.plot.config import *

def hrplot():

    physmerge = cksgaia.io.load_table(full_sample)
    crop = cksgaia.io.load_table(filtered_sample)

    fig = pl.figure(figsize=(12,8))

    pl.semilogy()

    allstars = physmerge.drop_duplicates(subset=['id_starname'])
    stars = crop.drop_duplicates(subset=['id_starname'])

    pl.plot(allstars['cks_steff'], allstars['giso_srad'], 'k.', alpha=0.2, markersize=8)
    pl.plot(stars['cks_steff'], stars['giso_srad'], 'bo', markersize=6)
    #pl.plot(gap['cks_steff'], gap['iso_srad'], 'ro', markersize=6)

    ls, li = 0.00025, 0.20
    xm = np.linspace(4000, 6750, 20)
    ym = 10**(ls*(xm-5500)+li)
    pl.plot(xm, ym, 'k--', lw=2)

    ax = pl.gca()
    ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%0.1f'))
    ax.yaxis.set_ticks(np.logspace(np.log10(0.5), np.log10(6.0), 8))
    #pl.yticks(np.logspace(np.log10(0.5), np.log10(6.0), 8))
    ax.minorticks_off()

    pl.xlim(6600, 4500)
    pl.ylim(0.5,6.0)

    pl.grid(lw=2)

    pl.xlabel('Stellar effective temperature [K]')
    pl.ylabel('Stellar radius [Solar radii]')


def filter_plot():
    matplotlib.rcParams['font.size'] = 18

    physmerge = cksgaia.io.load_table(full_sample)

    _ = cksgaia.io.apply_filters(physmerge, mkplot=True)


def simplehist(physmerge, color='k', annotate='err', nbins=36, fill_valley=True, aloc=(0.47, 0.85), eloc=(8.0, 85),
               stacked=False, weighted=False, va_anno=False, nstars=1, unc=False, clim=0.0, alpha=1):
    if weighted:
        weights = physmerge['weight'].values
    else:
        weights = 1 + np.zeros_like(physmerge['giso_prad'].dropna().values)

    xbins = np.logspace(np.log10(0.5), np.log10(20), nbins)
    N, edges = np.histogram(physmerge['giso_prad'].dropna().values, bins=xbins, weights=weights)
    Nd, edges = np.histogram(physmerge['giso_prad'].dropna().values, bins=xbins)
    centers = 0.5 * (edges[1:] + edges[:-1])
    p1loc = edges[(edges > 1.0) & (edges < 1.5)][np.argmax(N[(edges > 1.0) & (edges < 1.5)])]
    p2loc = edges[(edges > 2.0) & (edges < 4.0)][np.argmax(N[(edges > 2.0) & (edges < 4.0)])]
    p1 = np.max(N[(edges > 1.0) & (edges < 1.5)])
    p2 = np.max(N[(edges > 2.0) & (edges < 4.0)])
    top = np.min([p1, p2])
    bottom = np.min(N[(edges > 1.0) & (edges < 4.0)])
    vx = edges[(edges >= p1loc) & (edges <= p2loc)]
    vn = N[(edges >= p1loc) & (edges <= p2loc)]

    v = 0
    for i, (loc, val) in enumerate(zip(vx, vn)):
        if val < top:
            v += top - val
            try:
                if fill_valley:
                    pl.fill_between([loc, vx[i + 1]], top, val, color='0.5')
            except:
                pass

    v /= np.sum(N)
    v *= 100

    # radii limits for V_a calc
    v1, v2 = 1.64, 1.97
    valley = physmerge.query('giso_prad > %4.3f & giso_prad < %4.3f' % (v1, v2))
    # valley = np.where((physmerge['giso_prad'].values > v1) & (physmerge['giso_prad'].values < v2))[0]
    in_valley = len(valley)  # / (v2-v1)

    p1l, p1r = 1.2, 1.44
    peak1 = physmerge.query('giso_prad > %4.3f & giso_prad < %4.3f' % (p1l, p1r))
    # peak1 = np.where((physmerge['giso_prad'] > p1l).values & (physmerge['giso_prad'].values < p1r))[0]
    peak1 = len(peak1)  # / (p1r - p1l)

    p2l, p2r = 2.19, 2.62
    peak2 = physmerge.query('giso_prad > %4.3f & giso_prad < %4.3f' % (p2l, p2r))
    # peak2 = np.where((physmerge['giso_prad'].values > p2l) & (physmerge['giso_prad'].values < p2r))[0]
    peak2 = len(peak2)  # / (p2r - p2l)

    out_valley = np.mean([peak1, peak2])
    v = in_valley / out_valley
    # print peak1, in_valley, peak2, out_valley, v

    if va_anno:
        for l in [v1, v2]:
            pl.axvline(l, color='b', linestyle='dotted', lw=2)
        for l in [p1l, p1r] + [p2l, p2r]:
            pl.axvline(l, color='r', linestyle='dotted', lw=2)

    # print top, bottom, (top-bottom)/top, p1loc, p2loc, v

    pl.step(centers, N / nstars, where='mid', lw=3, color=color, alpha=alpha)
    bad = np.where(centers <= clim)[0]
    good = np.where(centers > clim)[0]
    if len(bad) > 0:
        pl.step(centers[bad], N[bad] / nstars, where='mid', color='0.8', lw=3)
    if unc:
        err = (np.sqrt(Nd) * (N / Nd)) / nstars
        err[np.isnan(err)] = 0
        _, caps, _ = pl.errorbar(centers[good], N[good] / nstars, fmt='.', yerr=err[good], lw=2, capsize=6, color=color)
        for cap in caps:
            cap.set_markeredgewidth(2)

    # physmerge['giso_prad'].hist(bins=xbins, histtype='step', lw=3, color=color)

    # pl.fill_between(edges[(edges >= p1loc) & (edges <= p2loc)],top, N[(edges >= p1loc) & (edges <= p2loc)],
    #                where=N[(edges >= p1loc) & (edges <= p2loc)] < top, color='0.5')

    x, y = eloc
    err1, err2 = cksgaia.misc.frac_err(physmerge, x, 'giso_prad')

    # if annotate == 'err':
    _, caps, _ = pl.errorbar([x], [y], fmt='.', ms=0.1, xerr=[[err1], [err2]], capsize=6, lw=2, color=color)
    for cap in caps:
        cap.set_markeredgewidth(2)
    pl.annotate("typical\nuncert.", xy=(x, y), xytext=(0, -50),
                xycoords="data", textcoords="offset points",
                horizontalalignment='center', fontsize=18)
    if annotate != 'err':
        pl.annotate(annotate, xy=aloc, xycoords='axes fraction', fontsize=20,
                    horizontalalignment='right')
        if fill_valley and va_anno:
            # pl.annotate("$V_A=%5.3f$" % v, xy=(aloc[0], aloc[1]-0.1), xycoords='axes fraction')
            pl.annotate("$V_A=%5.3f$" % v, xy=(0.50, aloc[1] - 0.11), xycoords='axes fraction')
    pl.xlim(0.7, 20)
    if not stacked:
        pl.ylabel('Number of Planets')
    pl.xlabel('Planet Size [Earth radii]')
    pl.semilogx()
    pl.xticks([0.7, 1.0, 1.5, 2.4, 3.5, 6, 10.0, 20.0])
    if stacked:
        ax = pl.gca()
        # yticks = ax.yaxis.get_major_ticks()
        # yticks[0].label1.set_visible(False)
        for label in ax.yaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        pl.xticks([])
    ax = pl.gca()
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%0.1f'))

    pl.grid(False)

    return v


def magcuts():
    physmerge = cksgaia.io.load_table(full_sample)

    figure = pl.figure(figsize=(10, 12))
    nrow = 3
    ncol = 1
    plti = 1

    xticks = [.7, 1.0, 1.3, 1.75, 2.4, 3.5, 4.5, 6]

    pl.subplot(nrow, ncol, plti)
    pl.subplots_adjust(hspace=0, top=0.98, bottom=0.10, left=0.19)

    bright = physmerge.query('kic_kepmag < 13.5')
    medium = physmerge.query('kic_kepmag < 14.2 & kic_kepmag >= 13.5')
    faint = physmerge.query('kic_kepmag >= 14.2')

    f1 = physmerge.query('kic_kepmag <= 14.2')
    f2 = physmerge.query('kic_kepmag <= 14.0')

    v = simplehist(bright, fill_valley=False, nbins=36, color='k', unc=True, aloc=(0.85, 0.80),
                                   annotate='$K_P < 13.5$', stacked=False, va_anno=False, eloc=(4.5, 50))

    ax = pl.gca()
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)
    for label in ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    pl.xticks([])

    pl.ylabel('')
    pl.ylim(0, 70)
    pl.xlim(0.7, 8)
    pl.xticks([])
    plti += 1

    pl.subplot(nrow, ncol, plti)
    v = simplehist(medium, fill_valley=False, nbins=36, color='k', unc=True, aloc=(0.85, 0.80),
                                   annotate='$13.5 < K_P \leq 14.2$', stacked=False,
                                   va_anno=False, weighted=False, eloc=(4.5, 45))

    ax = pl.gca()
    # yticks = ax.yaxis.get_major_ticks()
    # yticks[0].label1.set_visible(False)
    for label in ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    ax.yaxis.get_ticklabels()[0].set_visible(False)
    # ax.yaxis.get_ticklabels()[-1].set_visible(False)
    pl.xticks([])
    pl.ylabel('')
    pl.ylim(0, 70)
    pl.xlim(0.7, 8)
    pl.xticks([])
    plti += 1

    pl.ylabel('Number of Planets')

    pl.subplot(nrow, ncol, plti)
    v = simplehist(faint, fill_valley=False, nbins=36, color='k', unc=True, aloc=(0.85, 0.80),
                                   annotate='$K_P \geq 14.2$', stacked=False,
                                   va_anno=False, weighted=False, eloc=(4.5, 50))

    ax = pl.gca()
    # yticks = ax.yaxis.get_major_ticks()
    # yticks[0].label1.set_visible(False)
    for label in ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    pl.xticks([])

    pl.ylabel('')
    pl.ylim(0, 70)
    pl.xlim(0.7, 8)
    pl.xticks(xticks)


def depth_hist():
    physmerge = cksgaia.io.load_table(filtered_sample)

    fig = pl.figure(figsize=(10, 7))

    x = 7.0
    y = 75.0
    err1, err2 = cksgaia.misc.frac_err(physmerge, x, 'koi_ror')

    xbins = np.logspace(np.log10(0.1), np.log10(30), 40)
    pl.hist(physmerge['koi_ror'] * 100, bins=xbins, histtype='step', lw=4, color='k')
    # physmerge['koi_ror'].hist(bins=xbins, histtype='step', lw=4, color='k')
    _, caps, _ = pl.errorbar([x], [y], fmt='k.', xerr=[[err1], [err2]], capsize=6, lw=2, ms=0.1)
    for cap in caps:
        cap.set_markeredgewidth(2)
    pl.annotate("typical\nuncert.", xy=(x, y), xytext=(0, -60),
                xycoords="data", textcoords="offset points",
                horizontalalignment='center')
    pl.xlim(0.3, 30.0)
    pl.ylim(0, 100)
    pl.ylabel('Number of Planets')
    pl.xlabel('Planet-star radius ratio [%]')
    # pl.title('NEA')
    pl.semilogx()
    pl.xticks([0.3, 1.0, 3.0, 10.0, 30.0])
    ax = pl.gca()
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%0.1f'))

    pl.grid(False)


def srad_hist():
    physmerge = cksgaia.io.load_table(filtered_sample)

    fig = pl.figure(figsize=(10, 7))

    x = 2.2
    y = 100.0
    err1, err2 = cksgaia.misc.frac_err(physmerge, x, 'giso_srad')

    xbins = np.logspace(np.log10(0.6), np.log10(3.0), 20)
    physmerge['giso_srad'].hist(bins=xbins, histtype='step', lw=4, color='k')
    _, caps, _ = pl.errorbar([x], [y], fmt='k.', xerr=[[err1], [err2]], capsize=6, lw=2, ms=0.1)
    for cap in caps:
        cap.set_markeredgewidth(2)
    pl.annotate("typical\nuncert.", xy=(x, y), xytext=(0, -60),
                xycoords="data", textcoords="offset points",
                horizontalalignment='center')
    pl.xlim(0.6, 3.0)
    pl.ylim(0, 130)
    pl.ylabel('Number of Planets')
    pl.xlabel('Stellar Radius [Solar radii]')
    # pl.title('NEA')
    pl.semilogx()
    ax = pl.gca()
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%0.1f'))
    pl.xticks([0.6, 0.8, 1.0, 1.5, 2.0, 3.0])

    pl.grid(False)


def srad_err_hist():
    old = cksgaia.io.load_table('j17').groupby('id_starname').nth(0)
    new = cksgaia.io.load_table(full_sample).groupby('id_starname').nth(0)

    old['iso_srad_frac_err'] = (old['iso_srad_err1'] - old['iso_srad_err2'])/2. / old['iso_srad']
    new['giso_srad_frac_err'] = (new['giso_srad_err1'] - new['giso_srad_err2'])/2. / new['giso_srad']

    print len(old), len(new)

    fig = pl.figure(figsize=(6, 4))
    old['iso_srad_frac_err'] *= 100
    new['giso_srad_frac_err'] *= 100

    med_old = np.nanmedian(old['iso_srad_frac_err'])
    med_new = np.nanmedian(new['giso_srad_frac_err'])

    xbins = np.logspace(np.log10(0.5), np.log10(30.0), 30)
    old['iso_srad_frac_err'].hist(bins=xbins, histtype='step', lw=2, color='0.7')
    new['giso_srad_frac_err'].hist(bins=xbins, histtype='step', lw=2, color='k')

    pl.axvline(med_old, color='0.7', linestyle='dashed', lw=1)
    pl.axvline(med_new, color='k', linestyle='dashed', lw=1)

    pl.annotate("median = {:.0f}%".format(np.round(med_old)), xy=(med_old, 160), xycoords='data',
                xytext=(-15, 0), textcoords='offset points', rotation=90, verticalalignment='left')
    pl.annotate("median = {:.0f}%".format(np.round(med_new)), xy=(med_new, 120), xycoords='data',
                xytext=(-15, 0), textcoords='offset points', rotation=90, verticalalignment='left')

    pl.xlim(0.5, 30.0)
    # pl.ylim(0, 130)
    pl.ylabel('Number of Stars')
    pl.xlabel('Fractional Stellar Radius Uncertainty [%]')
    # pl.title('NEA')
    pl.semilogx()
    ax = pl.gca()
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%0.1f'))
    pl.xticks([0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 30.0])

    pl.grid(False)

    pl.legend(['Johnson+17', 'this work'], loc='upper left')


def prad_err_hist():
    old = cksgaia.io.load_table('j17')
    new = cksgaia.io.load_table(full_sample)

    fig = pl.figure(figsize=(6, 4))

    old['iso_prad_frac_err'] = cksgaia.errors.frac_err(old['iso_prad'], old['iso_prad_err1'], old['iso_prad_err2'])
    new['giso_prad_frac_err'] = cksgaia.errors.frac_err(new['giso_prad'], new['giso_prad_err1'], new['giso_prad_err2'])

    old['iso_prad_frac_err'] *= 100
    new['giso_prad_frac_err'] *= 100

    print len(old), len(new)

    med_old = np.nanmedian(old['iso_prad_frac_err'])
    med_new = np.nanmedian(new['giso_prad_frac_err'])

    xbins = np.logspace(np.log10(0.5), np.log10(50.0), 50)
    old['iso_prad_frac_err'].hist(bins=xbins, histtype='step', lw=2, color='0.7')
    new['giso_prad_frac_err'].hist(bins=xbins, histtype='step', lw=2, color='k')

    pl.axvline(med_old, color='0.7', linestyle='dashed', lw=2)
    pl.axvline(med_new, color='k', linestyle='dashed', lw=2)

    pl.annotate("median = {:.0f}%".format(np.round(med_old)), xy=(med_old, 200), xycoords='data',
                xytext=(-15,0), textcoords='offset points', rotation=90, verticalalignment='left')
    pl.annotate("median = {:.0f}%".format(np.round(med_new)), xy=(med_new, 200), xycoords='data',
                xytext=(-15,0), textcoords='offset points', rotation=90, verticalalignment='left')

    pl.xlim(0.5, 40.0)
    # pl.ylim(0, 130)
    pl.ylabel('Number of Planets')
    pl.xlabel('Fractional Planet Radius Uncertainty [%]')
    # pl.title('NEA')
    pl.semilogx()
    ax = pl.gca()
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%0.1f'))
    pl.xticks([1.0, 2.0, 3.0, 5.0, 10.0, 15.0, 30.0])

    pl.grid(False)

    pl.legend(['Johnson+17', 'this work'], loc='upper left')


def parallax_err_hist():
    old = cksgaia.io.load_table('iso-old')
    new = cksgaia.io.load_table(full_sample)

    fig = pl.figure(figsize=(12, 8))

    new['gaia2_sparallax_err'] *= 1e3

    med_new = np.nanmedian(new['gaia2_sparallax_err'])

    xbins = np.logspace(np.log10(1.0), np.log10(300.0), 80)
    new['gaia2_sparallax_err'].hist(bins=xbins, histtype='step', lw=4, color='k')

    pl.axvline(med_new, color='k', linestyle='dashed', lw=4)

    pl.annotate("median = {:.0f}$\mu$\"".format(np.round(med_new)), xy=(med_new, 40), xycoords='data',
                xytext=(-22, 0), textcoords='offset points', rotation=90, verticalalignment='left')

    pl.xlim(10.0, 300.0)
    # pl.ylim(0, 130)
    pl.ylabel('Number of Stars')
    pl.xlabel('Parallax Uncertainty [$\mu$ arcseconds]')
    # pl.title('NEA')
    pl.semilogx()
    ax = pl.gca()
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%0.1f'))
    pl.xticks([10.0, 30.0, 100.0, 300.0])

    pl.grid(False)

    # pl.legend(['Johnson+17', 'this work'], loc='upper left')