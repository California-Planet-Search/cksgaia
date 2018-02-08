import pylab as pl
import matplotlib
import numpy as np

import cksgaia.io


def hrplot():

    physmerge = cksgaia.io.load_table('cks-physical-merged')
    crop = cksgaia.io.load_table('cks-physical-filtered')

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
