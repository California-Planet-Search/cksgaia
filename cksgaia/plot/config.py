# put common plotting variables and values here
import matplotlib

matplotlib.rcParams['font.size'] = 24
matplotlib.rcParams['figure.figsize'] = (12,10)


full_sample = 'johnson17'
filtered_sample = 'fulton17'


print("Making plots with the {} table for the full sample and {} table for the filtered sample."
      .format(full_sample, filtered_sample))