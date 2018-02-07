import os
import cPickle as pickle

import pandas as pd
import numpy as np
from scipy.io import idl
from astropy.io import fits 
from astropy import constants as c 
from astropy import units as u

import cksmet.crossmatch
import cksmet.cuts
import cksmet.grid
import cksmet.analysis
import cksmet.pdplus
import cksspec.io
import cksmet.pdplus
import cksspec.utils

import cksmet.population


DATADIR = os.path.join(os.path.dirname(__file__),'../data/')
FLUX_EARTH = (c.L_sun / (4.0 * np.pi * c.au**2)).cgs.value
COMPLETENESS_FILE = os.path.join(DATADIR, 'comp.pkl')
from cksmet.calibrate import LAMO_SAMPLE, CAL_SAMPLE

def load_table(table, cache=0, cachefn='load_table_cache.hdf', verbose=False):
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

    elif table=='cksbin-fe':
        bins = [0.7, 1.0, 1.4, 2.0, 2.8, 4.0, 5.7, 8.0, 11.3, 16]
        cks = load_table('cks')
        #cks = cks.query(' -0.3 < feh_cks < 0.5')
        cksbin = table_bin(cks, bins)
        df = cksbin

    elif table=='cksbin-nplanets':
        cks = load_table('cks')
        g = cks.groupby('nplanets')
        g = g['feh_cks']
        dfbin = pd.DataFrame(index=g.first().index)
        dfbin['multiplicity'] = dfbin.index
        dfbin['nplanets'] = g.count() 
        dfbin['nstars'] = dfbin['nplanets']/ dfbin.multiplicity
        dfbin['fe_mean'] = g.mean()
        dfbin['fe_std'] = g.std()
        dfbin['fe_mean_err'] = dfbin['fe_std']/np.sqrt(dfbin['nstars'])
        df = dfbin

    elif table=='lamost':
        df = load_table('lamost-dr3')


    elif table=='lamost-dr3':
        dr3 = pd.read_hdf('data/lamost/dr3.hdf','dr3')
        stellar = load_table('mathur17',cache=1)

        dr3 = dr3.drop(['ra','dec'],axis=1)
        dr3 = dr3.rename(columns={'ra_obs':'ra','dec_obs':'dec'})
        s17 = cksmet.io.load_table('mathur17')
        s17 = s17['id_kic m17_ra m17_dec'.split()]
        s17 = s17.rename(columns={'m17_ra':'ra','m17_dec':'dec'})
        df =  cksmet.crossmatch.match_coords(dr3,s17,dist=1.2)
        df = df['id_kic teff teff_err logg logg_err feh feh_err snru snrg snrr snri snrz rv rv_err tsource tcomment tfrom'.split()]
        df['steff'] = df['teff']
        df['slogg'] = df['logg']
        df['smet'] = df['feh']
        df['steff_err1'] = df['teff_err']
        df['steff_err2'] = -1.0 * df['teff_err']
        df['slogg_err1'] = df['logg_err']
        df['slogg_err2'] = -1.0 * df['logg_err']
        df['slogg_err1'] = df['logg_err']
        df['slogg_err2'] = -1.0 * df['logg_err']
        df['smet_err1'] = df['feh_err']
        df['smet_err2'] = -1.0 * df['feh_err']
        df = add_prefix(df,'lamo_')
        df = pd.merge(df,stellar)

        # Select row with highest SNR
        df = df.sort_values(by=['id_kic','lamo_snrg'])
        df = df.groupby('id_kic',as_index=False).last()
        
        # Only work with stars brighter than 14 mag 
        df = df.query('m17_kepmag < 14.2')

    elif table=='lamost-cal':
        fn = os.path.join(DATADIR,'lamost-cal.hdf')
        df = pd.read_hdf(fn, 'lamost-cal')

    elif table=='cks':
        df = pd.read_csv('/Users/petigura/Research/CKS-Physical/data/cks_physical_merged.csv')

    elif table=='lamost+cks':
        cks = load_table('cks')
        lamost = load_table('lamost',cache=1)
        df = pd.merge(lamost,cks,on='id_kic')
        return df

    elif table=='cks-cuts':
        cks = load_table('cks')
        cks['id_koi'] = cks.id_starname.str.slice(start=1)
        cks['id_koi'] = pd.to_numeric(cks.id_koi,errors='coerce')

        m17 = load_table('mathur17')

        furlan2 = load_table('furlan-table2')
        furlan9 = load_table('furlan-table9')

        df = pd.merge(cks,m17)
        df = pd.merge(df,furlan2[['id_koi','furlan_ao_obs']],how='left')
        df = pd.merge(df,furlan9[['id_koi','furlan_rcorr_avg']],how='left')

        df = cksmet.cuts.add_cuts(df, cksmet.cuts.plnt_cuttypes, 'cks')

    elif table=='lamost-cks-calibration-sample-noclip':
        # Load up specmatch results specmatch sample
        cal = cksmet.io.load_table(CAL_SAMPLE,cache=1)
        cal = cksmet.io.sub_prefix(cal,'cks_')
        cal = cal['id_kic steff slogg smet svsini kic_kepmag'.split()]

        lamo = cksmet.io.load_table(LAMO_SAMPLE,cache=1)
        lamo = cksmet.io.sub_prefix(lamo,'lamo_')
        lamo = lamo['id_kic steff slogg smet snrg'.split()]
        namemap = {'steff':'teff','slogg':'logg','smet':'fe'}

        cal = cal.rename(columns=namemap)
        lamo = lamo.rename(columns=namemap)

        # Merge two catalogs and compute the differences
        df = cksspec.utils.merge_spec_tables(
            lamo, cal, on=['id_kic'], suffixes=['_new','_lib']
        )
        print ""
        print "_new =  {}".format(LAMO_SAMPLE)
        print "_lib =  {}".format(CAL_SAMPLE)
        print ""
        
        df = df.groupby('id_kic',as_index=False).first()

    elif table=='lamost-cks-calibration-sample':
        df = load_table('lamost-cks-calibration-sample-noclip')
        df = df.query('abs(fe_diff) < 0.2')
        print "{} stars in comparison after removing outliers".format(len(df))
        
    elif table=='lamost-cal-cuts':
        # Apply similar set of cuts to the lamost sample.
        df = load_table('lamost-cal')

        m17 = load_table('mathur17')
        df = pd.merge(df,m17)

        cuttypes = cksmet.cuts.lamo_cuttypes
        df = cksmet.cuts.add_cuts(df, cuttypes, 'lamo')


    elif table=='lamost-cal-cuts+cdpp':
        lamo = cksmet.io.load_table('lamost-cal-cuts',cache=1)
        huber14 = cksmet.io.load_table('huber14+cdpp',cache=1)
        huber14 = huber14['id_kic kepmag cdpp3 cdpp6 cdpp12'.split()]
        lamo = pd.merge(lamo,huber14,on='id_kic')
        df = lamo

    elif table=='field-cuts':
        df = load_table('mathur17+cdpp',cache=1)
        cuttypes = cksmet.cuts.field_cuttypes
        df = cksmet.cuts.add_cuts(df, cuttypes, 'field')

    elif table=='cdpp':
        fn = os.path.join(DATADIR,'kic_q0_q17.dat')
        df = idl.readsav(fn) 
        df = df['kic']
        df = cksmet.pdplus.LittleEndian(df) # Deals with the byte order
        df = pd.DataFrame(df)
        df = df.rename(columns={
            'KEPMAG':'kepmag','KICID':'id_kic','CDPP3':'cdpp3','CDPP6':'cdpp6',
            'CDPP12':'cdpp12'}
        )
        df = df['id_kic kepmag cdpp3 cdpp6 cdpp12'.split()]
        for col in 'cdpp3 cdpp6 cdpp12'.split():
            cdpp = np.vstack(df.ix[:,col])
            cdpp[cdpp==0.0] = np.nan
            cdpp = np.nanmedian(cdpp,axis=1)
            df[col] = cdpp
            df['log'+col] = np.log10(cdpp)
        
    elif table=='mathur17+cdpp':
        m17 = load_table('mathur17')
        cdpp = load_table('cdpp')
        df = pd.merge(m17,cdpp)

    elif table=='huber14+cdpp':
        m17 = load_table('mathur17')
        h14 = load_table('huber14')
        cdpp = load_table('cdpp')
        df = pd.merge(m17,cdpp)
        df = pd.merge(df,h14)

    elif table=='mathur17':
        df = cksspec.io.load_table('stellar17')
        namemap = {}
        for col in list(df.columns):
            if col[:3]=='kic':
                namemap[col] = col.replace('kic','m17')
        df = df.rename(columns=namemap)
        # add place holder column for quarter zero
        df['st_quarters'] = '0'+df.st_quarters 
        df['tobs'] = 0
        for q in range(18):
            # Was quarter observed
            qobs = df.st_quarters.str.slice(start=q,stop=q+1).astype(int)
            if q==17:
                qobs=0
            df['tobs'] += qobs * long_cadence_day * lc_per_quarter[q]

    elif table=='huber14':
        df = cksspec.io.load_table('huber14')
        

    elif table=='furlan-table2':
        tablefn = os.path.join(DATADIR,'furlan17/Table2.txt')
        df = pd.read_csv(tablefn,sep='\s+')
        namemap = {
            'KOI':'id_koi','KICID':'id_kic','Observatories':'ao_obs'
        }
        df = df.rename(columns=namemap)[namemap.values()]
        df['id_starname'] = ['K'+str(x).rjust(5, '0') for x in df.id_koi] # LMW convert id_koi to id_starname for merge with cks
        df = add_prefix(df,'furlan_')

    elif table=='furlan-table9':
        tablefn = os.path.join(DATADIR,'furlan17/Table9.txt')
        names = """
        id_koi hst hst_err i i_err 692 692_err lp600 lp600_err jmag jmag_err 
        kmag kmag_err jkdwarf jkdwarf_err jkgiant jkgiant_err rcorr_avg 
        rcorr_avg_err
        """.split()

        df = pd.read_csv(tablefn,sep='\s+',skiprows=2,names=names)
        df['id_starname'] = ['K'+str(x).rjust(5, '0') for x in df.id_koi] # LMW convert id_koi to id_starname for merge with cks
        df = add_prefix(df,'furlan_')

    elif table=='occur-surface':
        df = cksmet.surface.compute_occur_surface() 

    elif table=='per-prad-population':
        p, nplanets = cksmet.population.load_pinky() 
        df = cksmet.population.load_population(p, nplanets)

    elif table == 'cks+kmag+fakegaia':
        cks = load_table('cks+nea+iso-floor')
        cks['iso_sparallax_err1'] /= 5
        cks['iso_sparallax_err2'] /= 5
        cks['iso_sparallax_err1'] = np.sqrt(cks['iso_sparallax_err1']**2 + (3e-5)**2)  # 30 microarcsec floor
        cks['sio_sparallax_err2'] = -cks['iso_sparallax_err1']

        df = cks

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

def load_object(key, cache=1):
    pklfn = os.path.join(DATADIR,key+'.pkl')
    if cache==1:
        with open(pklfn,'r') as f:
            obj = pickle.load(f)
            return obj

    elif cache==2:
        obj = load_object(key,cache=0)
        with open(pklfn,'w') as f:
            pickle.dump(obj,f)
        
        return obj

    if key.count('occur')==1:
        obj = cksmet.analysis.load_occur(key)
    elif key.count('comp')==1:
        obj = cksmet.analysis.load_completeness()
    elif key.count('fit')==1:
        obj = cksmet.analysis.load_fit(key)
    return obj

def load_occur(key, cache=1):
    pklfn = os.path.join(DATADIR,key+'.pkl')
    if cache==1:
        with open(pklfn,'r') as f:
            fit = pickle.load(f)
            return fit 

    elif cache==2:
        fit = load_occur(key,cache=0)
        with open(pklfn,'w') as f:
            pickle.dump(fit,f)



def latex_table(table):
    tab = load_table(table)
    if table=='cksbin-fe':
        for i, row in tab.iterrows():
            line = "{bin0:.1f}--{bin1:.1f} &"
            line+=" {count:.0f} & "
            line+=" {fe_mean:+.3f} & "
            line+=" {fe_mean_err:.3f} &"
            line+=" {fe_p01:+.3f} &"
            line+=" {fe_p50:+.3f} &"
            line+=" {fe_p99:+.3f}"
            line+=r"\\"
            line = line.format(**row)
            print line

# Pulled from Kepler data characteristics
lc_per_quarter = {
    0:476,
    1:1639,
    2:4354,
    3:4370,
    4:4397,
    5:4634,
    6:4398,
    7:4375,
    8:3279,
    9:4768,
    10:4573,
    11:4754,
    12:4044,
    13:4421,
    14:4757,
    15:4780,
    16:4203,
    17:1556,
}
long_cadence_day = 29.7 / 60.0 / 24.0 # long cadence measurement in days



