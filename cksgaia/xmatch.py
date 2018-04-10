"""
Faciliate cross-matching between Gaia and the KIC
"""
import cksgaia.io
import pandas as pd
import numpy as np
max_sep = 1.0 # arcsec
def gaia1(df, gaia, key):
    """
    Crossmatch the sources in Gaia 1

    Args:
        df (pandas.DataFrame): Target catalog 
        gaia (pandas.DataFrame): Gaia DR1 table
        key (str): key to join on
    """

    assert len(df)==len(df.drop_duplicates()), "No duplicate stars"
    gaia['id_gaia1']= gaia.id_gaia.astype(str)  # handle nans

    # just want the stars
    m = pd.merge(df,gaia,on=key,how='left')
    m['id_gaia1']= m.id_gaia1.fillna('-99')
    m['id_gaia1']= m.id_gaia1.astype(np.int64)
    ndf = len(df)
    print "max(gaia1_angdist) = {} (arcsec)".format(m.gaia1_angdist.max())
    print "{} gaia sources within 8 arcsec of {} target sources".format(
        len(m),ndf
    )

    # count the number of stars within 8 arcsec
    m.index = m[key]
    g = m.groupby(key)
    m['gaia1_nsource'] = g.size()
    m['gaia1_gflux_sum'] = g['gaia1_gflux'].sum()
    m['gaia1_gflux_nearest'] = g['gaia1_gflux'].first()
    m['gaia1_gflux_ratio'] = m['gaia1_gflux_sum']/m['gaia1_gflux_nearest']

    m['gaia1_angdist_n1'] = g['gaia1_angdist'].nth(1) 
    g = m.groupby(key,as_index=False)
    m = g.first() # designate the closest source as the correct source


    #ratio = 1.05
    #print "{}/{} sources have gaia dilution of < {}".format(
    #    (m.gaia1_gflux_ratio < ratio).sum(),ndf,ratio
    #)
    return m
