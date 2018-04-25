# put common plotting variables and values here
import matplotlib
import pylab as pl
import os

import cksgaia

matplotlib.rcParams['font.size'] = 24
matplotlib.rcParams['figure.figsize'] = (12, 10)

afs = 24  # font size for annotations

full_sample = 'j17+m17+gaia2'
filtered_sample = 'fulton17-weights'

modpath = '/'.join(os.path.dirname(cksgaia.__file__).split('/')[:-1])

print("Making plots with the {} table for the full sample and {} table for the filtered sample."
      .format(full_sample, filtered_sample))


def logscale_rad():
      pl.semilogx()

      pl.xlabel('Planet Size [Earth radii]')

      pl.xticks([.7, 1.0, 1.3, 1.8, 2.4, 3.5, 4.5, 6, 8.0])
      # pl.xticks(np.logspace(np.log10(0.7), np.log10(8), 6))

      ax = pl.gca()
      ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
      ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%0.1f'))

      pl.xlim(0.7, 6.0)

      return ax
