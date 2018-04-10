# put common plotting variables and values here
import matplotlib
import os

import cksgaia

matplotlib.rcParams['font.size'] = 24
matplotlib.rcParams['figure.figsize'] = (12, 10)

afs = 24  # font size for annotations

full_sample = 'j17'
filtered_sample = 'fulton17-weights'

modpath = '/'.join(os.path.dirname(cksgaia.__file__).split('/')[:-1])

print("Making plots with the {} table for the full sample and {} table for the filtered sample."
      .format(full_sample, filtered_sample))