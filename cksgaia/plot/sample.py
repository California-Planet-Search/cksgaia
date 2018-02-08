import pylab as pl
import matplotlib
import numpy as np

import cksgaia.io
import cksgaia.misc

def hrplot():

    physmerge = cksgaia.io.load_table('johnson17')
    crop = cksgaia.io.load_table('fulton17')

    fig = pl.figure(figsize=(12,8))

    pl.semilogy()

    allstars = physmerge.drop_duplicates(subset=['id_starname'])
    stars = crop.drop_duplicates(subset=['id_starname'])

    pl.plot(allstars['cks_steff'], allstars['iso_srad'], 'k.', alpha=0.2, markersize=8)
    pl.plot(stars['cks_steff'], stars['iso_srad'], 'bo', markersize=6)
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
    pl.minorticks_off()

    pl.xlim(6600, 4500)
    pl.ylim(0.5,6.0)

    pl.grid(lw=2)

    pl.xlabel('Stellar effective temperature [K]')
    pl.ylabel('Stellar radius [Solar radii]')


def filter_plot():
    matplotlib.rcParams['font.size'] = 18

    physmerge = cksgaia.io.load_table('johnson17')

    _ = cksgaia.io.apply_filters(physmerge, mkplot=True)


def simplehist(physmerge, color='k', annotate='err', nbins=36, fill_valley=True, aloc=(0.47, 0.85), eloc=(8.0, 85),
               stacked=False, weighted=False, va_anno=False, nstars=1, unc=False, clim=0.0, alpha=1):
    if weighted:
        weights = physmerge['weight'].values
    else:
        weights = 1 + np.zeros_like(physmerge['iso_prad'].dropna().values)

    xbins = np.logspace(np.log10(0.5), np.log10(20), nbins)
    N, edges = np.histogram(physmerge['iso_prad'].dropna().values, bins=xbins, weights=weights)
    Nd, edges = np.histogram(physmerge['iso_prad'].dropna().values, bins=xbins)
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
    valley = physmerge.query('iso_prad > %4.3f & iso_prad < %4.3f' % (v1, v2))
    # valley = np.where((physmerge['iso_prad'].values > v1) & (physmerge['iso_prad'].values < v2))[0]
    in_valley = len(valley)  # / (v2-v1)

    p1l, p1r = 1.2, 1.44
    peak1 = physmerge.query('iso_prad > %4.3f & iso_prad < %4.3f' % (p1l, p1r))
    # peak1 = np.where((physmerge['iso_prad'] > p1l).values & (physmerge['iso_prad'].values < p1r))[0]
    peak1 = len(peak1)  # / (p1r - p1l)

    p2l, p2r = 2.19, 2.62
    peak2 = physmerge.query('iso_prad > %4.3f & iso_prad < %4.3f' % (p2l, p2r))
    # peak2 = np.where((physmerge['iso_prad'].values > p2l) & (physmerge['iso_prad'].values < p2r))[0]
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

    # physmerge['iso_prad'].hist(bins=xbins, histtype='step', lw=3, color=color)

    # pl.fill_between(edges[(edges >= p1loc) & (edges <= p2loc)],top, N[(edges >= p1loc) & (edges <= p2loc)],
    #                where=N[(edges >= p1loc) & (edges <= p2loc)] < top, color='0.5')

    x, y = eloc
    err1, err2 = cksgaia.misc.frac_err(physmerge, x, 'iso_prad')

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
