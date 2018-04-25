import os

DATADIR = os.path.join(os.path.dirname(__file__),'../data/')

# Path to where the isochrones are stored
ISO_CSVFN = os.path.join(DATADIR,'isochrones_gaia2.csv')
CKS_CSVFN = os.path.join(DATADIR,'cks_v2.csv')
# MERGED_TABLE = os.path.join(DATADIR,'cks_fakegaia_merged.csv')
MERGED_TABLE = os.path.join(DATADIR,'cks_iso_gaia2_merged.csv')
MERGED_TABLE_NAME = os.path.join('cksgaia')

# slope and intercept for subgiant filter
ls, li = 0.00025, 0.20

# Total number of stars in the search sample
num_stars = 36075
