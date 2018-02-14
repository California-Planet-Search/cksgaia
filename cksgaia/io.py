import os

import pandas as pd
import numpy as np
import pylab as pl
import ebf

import cksgaia.plot
import cksgaia.completeness
from cksgaia.config import *

DATADIR = os.path.join(os.path.dirname(__file__), '../data/')

def load_table(table, cache=1, cachefn='load_table_cache.hdf', verbose=False):
    """Load tables used in cksmet

    Args:
        table (str): name of table. must be one of
            - nea 


        cache (Optional[int]): whether or not to use the cache
            - 0: don't use the cache recreate all files
            - 1: read from cache
            - 2: write tables to cache

    Returns:
        pandas.DataFrame: table

    """
    if cache==1:
        try:
            df = pd.read_hdf(cachefn,table)
            print "read table {} from {}".format(table,cachefn)
            return df
        except IOError:
            print "Could not find cache file: %s" % cachefn
            print "Building cache..."
            cache=2
        except KeyError:
            print "Cache not built for table: %s" % table
            print "Building cache..."
            cache=2

    if cache==2:
        df = load_table(table, cache=False)
        print "writing table {} to cache".format(table)
        df.to_hdf(cachefn,table)
        return df

    elif table=='stellar17':
        tablefn = os.path.join(DATADIR, 'kepler_stellar17.csv.gz')
        df = pd.read_csv(tablefn,sep='|',dtype={'st_quarters':'str'})
        namemap = {
            'kepid':'id_kic','kepmag':'kic_kepmag', 'teff': 'kic_steff',
            'st_quarters':'st_quarters','mass':'kic_smass',
            'st_radius':'kic_srad', 'jmag':'kic_jmag',
            'jmag_err':'kic_jmag_err','hmag':'kic_hmag',
            'hmag_err':'kic_hmag_err','kmag':'kic_kmag',
            'kmag_err':'kic_kmag_err'
        }
        df = df.rename(columns=namemap)[namemap.values()]

    elif table=='mathur17':
        df = cksgaia.io.load_table('stellar17')
        namemap = {}
        for col in list(df.columns):
            if col[:3]=='kic':
                namemap[col] = col.replace('kic','m17')
        df = df.rename(columns=namemap)

    elif table=='johnson17':
        df = pd.read_csv('data/cks_physical_merged.csv',index_col=0)

    elif table=='fulton17':
        df = load_table('johnson17')
        df = apply_filters(df)

    elif table=='fulton17-weights':
        df = load_table('johnson17')
        df = apply_filters(df)
        df = cksgaia.completeness.weight_merge(df)

    elif table=='cks-physical-merged+mathur17':
        df = load_table('cks-physical-merged')

    elif table=='j17+m17':
        df = load_table('johnson17')
        m17 = load_table('mathur17')
        df = pd.merge(df, m17, on='id_kic')

    elif table == 'j17+m17-fakegaia':
        df = load_table('j17+m17')
        iso_err = df.iso_sparallax_err1 / 5 
        gaia_err = np.sqrt(iso_err**2 + 0.02**2) # 20uas error floor
        df['iso_sparallax_err1'] = gaia_err
        df['iso_sparallax_err2'] = -1.0 * gaia_err

    elif table == 'kic':
        fname = os.path.join(DATADIR, 'kic_q0_q17.hdf')
        kic = pd.read_hdf(fname)
        s17 = load_table('stellar17')
        kicmerge = pd.merge(s17, kic, left_on='id_kic', right_on='KICID')
        df = kicmerge

    elif table == 'kic-filtered':
        df = load_table('kic')

        kicselect = df.query('kic_kepmag <= 14.2 & kic_steff >= 4700 & kic_steff <= 6500')
        kicselect = kicselect[kicselect['kic_srad'] <= 10**(ls*(kicselect['kic_steff']-5500)+li)]
        kicselect = kicselect.dropna(subset=['kic_smass']).reset_index()

        df = kicselect

    else:
        assert False, "table {} not valid table name".format(table)
    return df

def add_prefix(df,prefix,ignore=['id']):
    namemap = {}
    for col in list(df.columns):
        skip=False
        for _ignore in ignore:
            if col.count(_ignore) > 0:
                skip = True
        if not skip:
            namemap[col] = prefix + col 
    df = df.rename(columns=namemap)
    return df


def sub_prefix(df, prefix,ignore=['id']):
    namemap = {}
    for col in list(df.columns):
        skip=False
        for _ignore in ignore:
            if col.count(_ignore) > 0:
                skip = True
        if not skip:
            namemap[col] = col.replace(prefix,'') 
    df = df.rename(columns=namemap)
    return df

def load_mist():
    model = ebf.read(os.path.join(DATADIR,'mesa.ebf'))
    return model


def apply_filters(physmerge, mkplot=False, verbose=False, textable=False):
    # slope and intercept for subgiant filter
    ls, li = 0.00025, 0.20

    def _bipanel(physmerge, nrow, ncol, plti, eloc=(12.0, 100), aloc=(0.95, 0.85), atxt='full CKS sample',
                 letters=True):

        if mkplot:
            left = plti
            right = plti + 1

            letters = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
            left_letter = letters[left - 1]
            right_letter = letters[right - 1]

            nea = physmerge.copy()
            nea['iso_prad'] = nea['koi_prad']
            nea['iso_prad_err1'] = nea['koi_prad_err1']

            pl.subplot(nrow, ncol, left)
            v_cks = cksgaia.plot.sample.simplehist(physmerge, annotate='%s (%d)' % (atxt, len(physmerge)),
                                               stacked=True, va_anno=False,
                                               eloc=eloc, aloc=aloc, fill_valley=False)
            if letters:
                pl.annotate(left_letter + ")", xy=(0.03, 0.9),
                            xycoords='axes fraction', weight='bold')
            pl.xlabel('')
            ymax_left = pl.ylim()[1]
            ax = pl.gca()
            plti += 1

            texline = "%s  &  %4.3f  \\\\" % (atxt, v_cks)
            if textable:
                f = open('tmp.tex', 'a')
                print >>f, texline
                f.close()

            print texline

            ylimits = [0, ymax_left]
            pl.subplot(nrow, ncol, left)
            pl.ylim(ylimits)

            return left + 1

    nrow = 7
    ncol = 1
    plti = 1

    if mkplot:
        figure = pl.figure(figsize=(8, 20))
        pl.subplots_adjust(hspace=0, wspace=0.25, top=0.98, bottom=0.05, right=0.96, left=0.15)

    if verbose:
        print "Initial catalogue = %d planets." % len(physmerge)

    plti = _bipanel(physmerge, nrow, ncol, plti)

    pre = len(physmerge)
    crop = physmerge[physmerge['cks_fp'] == False]
    post = len(crop)
    if verbose:
        print "False positive filter removes %d planets." % (pre - post)
    plti = _bipanel(crop, nrow, ncol, plti, atxt='false pos. removed')

    pre = len(crop)
    crop = crop[crop['kic_kepmag'] <= 14.2]
    post = len(crop)
    if verbose:
        print "Kp < 14.2 filter removes %d planets." % (pre - post)
    plti = _bipanel(crop, nrow, ncol, plti, eloc=(12.0, 70), atxt='$Kp < 14.2$')

    # pre = len(crop)
    # crop = crop[crop['furlan_rcorr_avg'] <= 1.05]
    # if mkplot:
    #     pl.subplot(nrow,ncol,plti)
    #     v = cksrad.plotting.simplehist(crop, annotate='dilution $\leq$ 5%% (%d)' % len(crop), stacked=True, va_anno=True)
    # plti += 1
    # post = len(crop)
    # if verbose:
    #     print "Furlan+17 Rp correction < 5%% filter removes %d planets." % (pre-post)

    pre = len(crop)
    crop = crop[crop['koi_impact'] <= 0.7]
    post = len(crop)
    if verbose:
        print "b < 0.7 filter removes %d planets." % (pre - post)
    plti = _bipanel(crop, nrow, ncol, plti, eloc=(12.0, 60), atxt='$b < 0.7$')

    # pre = len(crop)
    # crop = crop[crop['vshaped'] >= 0.3]
    # pl.subplot(nrow,ncol,plti)
    # v = cksrad.plotting.simplehist(crop, annotate='v-shaped $\geq 0.3$ (%d)' % len(crop))
    # plti += 1
    # post = len(crop)
    # print "v-shaped >= 0.3 filter removes %d planets." % (pre-post)

    # pre = len(crop)
    # crop = crop[crop['iso_prad_err1']/crop['iso_prad'] <= 0.12]
    # if mkplot:
    #     pl.subplot(nrow,ncol,plti)
    #     v = cksrad.plotting.simplehist(crop, annotate='$\sigma(R_p)/R_p \leq 12$%% (%d)' % len(crop), stacked=True, va_anno=True)
    # plti += 1
    # post = len(crop)
    # if verbose:
    #     print "Rp_err < 12%% filter removes %d planets." % (pre-post)

    pre = len(crop)
    crop = crop[(crop['koi_period'] < 100)]  # & (crop['koi_period'] > 1)]
    post = len(crop)
    if verbose:
        print "Period filter removes %d planets." % (pre - post)
    plti = _bipanel(crop, nrow, ncol, plti, eloc=(12.0, 60), atxt='$P < 100$ d')

    pre = len(crop)
    crop = crop[crop['iso_srad'] <= 10 ** (ls * (crop['cks_steff'] - 5500) + li)]
    post = len(crop)
    if verbose:
        print "Subgiant filter removes %d planets." % (pre - post)
    plti = _bipanel(crop, nrow, ncol, plti, eloc=(12.0, 60), atxt='giant stars removed')

    pre = len(crop)
    crop = crop[(crop['cks_steff'] <= 6500) & (crop['cks_steff'] >= 4700)]
    post = len(crop)
    if verbose:
        print "Teff filter removes %d planets." % (pre - post)
    plti = _bipanel(crop, nrow, ncol, plti, eloc=(12.0, 60), aloc=(0.95, 0.85), atxt='4700 K < $T_{\\rm eff}$ < 6500 K')

    # pre = len(crop)
    # crop = crop[(crop['koi_snr'] >= 12)]
    # post = len(crop)
    # if verbose:
    #     print "SNR > 10 filter removes %d planets." % (pre-post)
    # plti = _bipanel(crop, nrow, ncol, plti, eloc=(8.0, 50), atxt='SNR > 10')

    if mkplot:
        # pl.subplot(nrow, ncol, plti-2)
        # pl.xticks([0.7, 1.0, 1.8, 3.5, 6.0, 10.0, 20.0])
        # pl.xlabel('Planet Size [Earth radii]', fontsize=28)

        pl.subplot(nrow, ncol, plti - 1)
        pl.xticks([0.7, 1.0, 1.8, 3.5, 6.0, 10.0, 20.0])
        pl.xlabel('Planet Size [Earth radii]', fontsize=28)

        # pl.annotate('Planet Size [Earth radii]', xy=(0.5, 0.03), xycoords='figure fraction',
        #                fontsize=28, horizontalalignment='center')
        pl.annotate('Number of Planets', xy=(0.04, 0.6), xycoords='figure fraction',
                    fontsize=28, horizontalalignment='center', rotation=90)
        # pl.annotate('Number of Planets', xy=(0.53, 0.6), xycoords='figure fraction',
        #                  fontsize=28, horizontalalignment='center', rotation=90)

    if verbose:
        print
        print "Final sample = %d planets." % len(crop)



    return crop
