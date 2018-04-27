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
        df = df.rename(columns=namemap)

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
        df = df.rename(columns=namemap)

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
        df = df.rename(columns=namemap)


    if mode.count('archive')==1:
        df['angdist'] *= 60*60

    df = df[namemap.values()]
    df = cksgaia.io.add_prefix(df, gaiadr+'_')
    return df

def read_xmatch_gaia2(fn):
    df = pd.read_csv(fn)
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
    df['steff'] = df['teff_val']
    df['steff_err1'] = df.eval('teff_percentile_upper - teff_val')
    df['steff_err2'] = df.eval('teff_percentile_lower - teff_val')
    df['srad'] = df['radius_val']
    df['srad_err1'] = df.eval('radius_percentile_upper - radius_val')
    df['srad_err2'] = df.eval('radius_percentile_lower - radius_val')
    df = df.rename(columns=namemap)
    df['angdist'] *= 60*60
    cols = namemap.values() + 'steff steff_err1 steff_err2 srad srad_err1 srad_err2 '.split()
    df = df[cols]
    df = cksgaia.io.add_prefix(df, 'gaia2_')
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
    Crossmatch the sources in Gaia 2

    Args:
        df (pandas.DataFrame): Target catalog 
        gaia (pandas.DataFrame): Gaia DR2 table
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
    print "max(gaia_angdist) = {} (arcsec)".format(m[gaiadr+'_angdist'].max())
    print "{} gaia sources within 8 arcsec of {} target sources".format(
        len(m),ndf
    )

    # count the number of stars within 8 arcsec
    m.index = m.id_kic
    g = m.groupby('id_kic')
    m[gaiadr+'_gflux_sum'] = g[gaiadr+'_gflux'].sum()
    m['absdiff_gmag_kepmag'] = np.abs(m['gaia2_gmag'] - m['kic_kepmag'])
    m['gaia2_n_8arcsec'] = g.size()

    # Match candidate
    mbest = m.query('gaia2_angdist < 1 and abs(kic_kepmag - gaia2_gmag) < 0.5') # within 1 arcsec
    mbest['gaia2_n_1arcsec'] = mbest.groupby('id_kic').size()

    print "{} gaia sources within 1 arcsec of {} target sources".format(
        len(mbest),mbest.id_kic.drop_duplicates().count()
    )

    mbest = mbest.sort_values(by=['id_kic','absdiff_gmag_kepmag'])
    g = mbest.groupby('id_kic',as_index=False)
    mbest['gaia2_n_1arcsec'] = g.size()
    mbest = g.nth(0) 
    mbest['gaia2_gflux_ratio'] = mbest.eval('gaia2_gflux_sum / gaia2_gflux')
    return mbest, m 
