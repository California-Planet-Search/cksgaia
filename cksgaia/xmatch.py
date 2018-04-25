"""
Faciliate cross-matching between Gaia and the KIC
"""
import cksgaia.io
import pandas as pd
import numpy as np
max_sep = 1.0 # arcsec

import astropy.io
import astropy.table
import os
from cksgaia.io import DATADIR

def create_xmatch_table():
    """
    Create VOTable and CSV tables of KIC and CKS samples
    """

    # CKS
    df = cksgaia.io.load_table('j17+m17')
    df = df['id_kic m17_ra m17_dec'.split()]
    df = df.groupby('id_kic',as_index=False).first()
    basename = "cks-{}".format(len(df))
    write_xmatch_table(df,basename)

    # KIC
    df = cksgaia.io.load_table('m17')
    df = df['id_kic m17_ra m17_dec'.split()]
    df = df.groupby('id_kic',as_index=False).first()
    basename = "m17-{}".format(len(df))
    write_xmatch_table(df,basename)

def write_xmatch_table(df,basename):
    """
    Create VOTable and CSV tables of KIC and CKS samples
    """
    fn = os.path.join(DATADIR,basename+'.csv')
    df.to_csv(fn,index=False)
    print "created {}".format(fn)

    fn = fn.replace('csv','vot')
    tab = astropy.table.Table.from_pandas(df)
    vot = astropy.io.votable.from_table(tab)
    astropy.io.votable.writeto(vot, fn)
    print "created {}".format(fn)

def read_xmatch_results(fn, mode):
    gaiadr, source = mode.split('-')
    df = pd.read_csv(fn)
    if mode=='gaia1-archive':
        namemap = {
            'source_id':'id_gaia1',
            'ra':'ra', 
            'dec':'dec',
            'parallax':'sparallax', 
            'parallax_error':'sparallax_err', 
            'phot_g_mean_flux':'gflux',
            'phot_g_mean_flux_error':'gflux_err',
            'phot_g_mean_mag':'gmag',
            'id_kic':'id_kic',
            'dist':'angdist',
        }

    elif mode=='gaia1-xmatch':
        namemap = {
            'source_id':'id_gaia1',
            'ra_ep2000':'ra', 
            'dec_ep2000':'dec',
            'parallax':'sparallax', 
            'parallax_error':'sparallax_err', 
            'phot_g_mean_flux':'gflux',
            'phot_g_mean_flux_error':'gflux_err',
            'phot_g_mean_mag':'gmag',
            'id_kic':'id_kic',
            'angDist':'angdist'
        }

    elif mode=='gaia2-archive':
        namemap = {
            'source_id':'id_gaia2',
            'ra':'ra', 
            'dec':'dec',
            'parallax':'sparallax', 
            'parallax_error':'sparallax_err', 
            'phot_g_mean_flux':'gflux',
            'phot_g_mean_flux_error':'gflux_err',
            'phot_g_mean_mag':'gmag',
            'parallax_over_error':'sparallax_over_err',
            'id_kic':'id_kic',
            'dist':'angdist',
        }
    elif gaiadr=='gaia2':
        pass

    df = df.rename(columns=namemap)
    if mode.count('archive')==1:
        df['angdist'] *= 60*60

    df = df[namemap.values()]
    df = cksgaia.io.add_prefix(df, gaiadr+'_')
    return df

def xmatch_gaia(df, gaia, key, gaiadr):
    """
    Crossmatch the sources in Gaia 1

    Args:
        df (pandas.DataFrame): Target catalog 
        gaia (pandas.DataFrame): Gaia DR1 table
        key (str): key to join on
        gaiadr (str): {'gaiadr1','gaiadr2'}
    """

    id_gaia = 'id_{}'.format(gaiadr)
    assert len(df.id_kic)==len(df.id_kic.drop_duplicates()), "No duplicate stars"
    gaia[id_gaia]= gaia[id_gaia].astype(str)  # handle nans

    # just want the stars
    m = pd.merge(df,gaia,on=key,how='left')
    m[id_gaia]= m[id_gaia].fillna('-99')
    m[id_gaia]= m[id_gaia].astype(np.int64)
    ndf = len(df)
    print "max(gaia1_angdist) = {} (arcsec)".format(m[gaiadr+'_angdist'].max())
    print "{} gaia sources within 8 arcsec of {} target sources".format(
        len(m),ndf
    )

    # count the number of stars within 8 arcsec
    m.index = m[key]
    g = m.groupby(key)
    m[gaiadr+'_nsource'] = g.size()
    m[gaiadr+'_gflux_sum'] = g[gaiadr+'_gflux'].sum()
    m[gaiadr+'_gflux_nearest'] = g[gaiadr+'_gflux'].first()
    m[gaiadr+'_gflux_ratio'] = m[gaiadr+'_gflux_sum']/m[gaiadr+'_gflux_nearest']
    m[gaiadr+'_angdist_n1'] = g[gaiadr+'_angdist'].nth(1) 
    m = m.sort_values(by=[key,gaiadr+'_angdist'])
    g = m.groupby(key,as_index=False)
    m = g.nth(0) # first is slow when there is a string in the columns
    return m


def xmatch_gaia2(df, gaia, key, gaiadr):
    """
    Crossmatch the sources in Gaia 1

    Args:
        df (pandas.DataFrame): Target catalog 
        gaia (pandas.DataFrame): Gaia DR1 table
        key (str): key to join on
        gaiadr (str): {'gaiadr1','gaiadr2'}
    """

    id_gaia = 'id_{}'.format(gaiadr)
    assert len(df.id_kic)==len(df.id_kic.drop_duplicates()), "No duplicate stars"
    gaia[id_gaia]= gaia[id_gaia].astype(str)  # handle nans

    # just want the stars
    m = pd.merge(df,gaia,on=key,how='left')
    m[id_gaia]= m[id_gaia].fillna('-99')
    m[id_gaia]= m[id_gaia].astype(np.int64)
    ndf = len(df)
    print "max(gaia1_angdist) = {} (arcsec)".format(m[gaiadr+'_angdist'].max())
    print "{} gaia sources within 8 arcsec of {} target sources".format(
        len(m),ndf
    )

    # count the number of stars within 8 arcsec
    m.index = m.id_gaia2
    g = m.groupby('id_kic')
    m['id_gaia2_best'] = False
    m['absdiff_gmag_kepmag'] = m['gaia2_gmag'] - m['kic_kepmag']
    for id_kic, idx in g.groups.iteritems():
        stars = m.loc[idx]
        stars.sort_values()
        if len(stars) > 4:
            
            import pdb;pdb.set_trace()
            
            
        print len(stars)
        

    #m[gaiadr+'_nsource'] = g.size()
    #m[gaiadr+'_gflux_sum'] = g[gaiadr+'_gflux'].sum()
    #m[gaiadr+'_gflux_nearest'] = g[gaiadr+'_gflux'].first()
    #m[gaiadr+'_gflux_ratio'] = m[gaiadr+'_gflux_sum']/m[gaiadr+'_gflux_nearest']
    #m[gaiadr+'_angdist_n0'] = g[gaiadr+'_angdist'].nth(0) 
    #m[gaiadr+'_angdist_n1'] = g[gaiadr+'_angdist'].nth(1) 
    #m = m.sort_values(by=[key,gaiadr+'_angdist'])
    #g = m.groupby(key,as_index=False)
    #m = g.nth(0) # first is slow when there is a string in the columns
    return m
