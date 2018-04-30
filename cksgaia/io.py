from collections import OrderedDict
import re

import pandas as pd
import numpy as np
import ebf
import cksgaia.plot
import cksgaia.completeness
from cksgaia.config import *
import cksgaia.extinction
import cksgaia.xmatch
import cksgaia.calc
from cksgaia.sample import apply_filters
from astropy import units as u
import astropy.io.ascii
DATADIR = os.path.join(os.path.dirname(__file__), '../data/')

def load_table(table, cache=1, cachefn='load_table_cache.hdf', verbose=False):
    """Load tables used in cksgaia

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

    if table=='coldefs':
        tablefn = os.path.join(DATADIR,'column-definitions.txt')
        colspecs = [(0,1),(3,4)]
        df = pd.read_fwf(
            tablefn, comment='#', widths=[20,100],
            names=['column','description']
        )

    # Mathur 2017
    elif table=='stellar17':
        tablefn = os.path.join(DATADIR, 'kepler_stellar17.csv.gz')
        df = pd.read_csv(tablefn,sep='|',dtype={'st_quarters':'str'})
        namemap = {
            'kepid':'id_kic','kepmag':'kic_kepmag', 'teff': 'kic_steff',
            'st_quarters':'st_quarters','mass':'kic_smass',
            'st_radius':'kic_srad', 'jmag':'kic_jmag',
            'jmag_err':'kic_jmag_err','hmag':'kic_hmag',
            'hmag_err':'kic_hmag_err','kmag':'kic_kmag',
            'kmag_err':'kic_kmag_err',
            'degree_ra':'kic_ra', 'degree_dec':'kic_dec'
        }
        df = df.rename(columns=namemap)[namemap.values()]

    elif table=='m17':
        df = load_table('stellar17')
        namemap = {}
        for col in list(df.columns):
            if col[:3]=='kic':
                namemap[col] = col.replace('kic','m17')
        df = df.rename(columns=namemap)

    # Gaia DR2
    elif table=='gaia2':
        fn = os.path.join(DATADIR, 'xmatch_m17_gaiadr2-result.csv')
        df = cksgaia.xmatch.read_xmatch_gaia2(fn)

    # Johnson 2017
    elif table=='j17':
        fn = MERGED_TABLE_OLD
        df = pd.read_csv(fn, index_col=0)

    # Silva 2015
    elif table=='silva15':
        fn = os.path.join(DATADIR,'silva15/silva-aguirre15.tex')
        df = read_silva15(fn)

    # Huber 2013
    elif table=='huber13':
        fn = os.path.join(DATADIR,'huber13/J_ApJ_767_127/table2.dat')
        readme = os.path.join(DATADIR,'huber13/J_ApJ_767_127/ReadMe')
        df = read_huber13(fn,readme)

    # Furlan 2017
    elif table=='furlan17-table2':
        fn = os.path.join(DATADIR,'furlan17/Table2.txt')
        df = read_furlan17_table2(fn)

    elif table=='furlan17-table9':
        fn = os.path.join(DATADIR,'furlan17/Table9.txt')
        df = read_furlan17_table9(fn)

    elif table=='fur17':
        tab2 = load_table('furlan17-table2')
        tab9 = load_table('furlan17-table9')
        cols = 'id_koi fur17_rcorr_avg fur17_rcorr_avg_err'.split()
        df = pd.merge(tab2,tab9[cols],how='left',on='id_koi')


    # Merged tables
    elif table=='m17+gaia2':
        print "performing crossmatch on gaia2"
        df = cksgaia.io.load_table('m17')
        df = df.rename(columns={'m17_kepmag':'kic_kepmag'})
        gaia = load_table('gaia2')
        stars = df['id_kic kic_kepmag'.split()].drop_duplicates()
        mbest,mfull = cksgaia.xmatch.xmatch_gaia2(stars,gaia,'id_kic','gaia2')
        df = pd.merge(df,mbest.drop(['kic_kepmag'],axis=1),on='id_kic')

    elif table=='m17+gaia2+j17':
        print "performing crossmatch on gaia2"
        df1 = cksgaia.io.load_table('m17+gaia2') 
        df2 = cksgaia.io.load_table('j17') 
        df2 = df2.drop(['kic_kepmag'],axis=1) # duplicated
        df = pd.merge(df1,df2,on='id_kic')

    elif table == 'kic':
        fname = os.path.join(DATADIR, 'kic_q0_q17.hdf')
        kic = pd.read_hdf(fname)
        s17 = load_table('stellar17')
        df = pd.merge(s17, kic, left_on='id_kic', right_on='KICID')

    elif table == 'kic-filtered':
        df = load_table('kic')
        df = df.query('kic_kepmag <= 14.2 & kic_steff >= 4700 & kic_steff <= 6500')
        df = df[df['kic_srad'] <= 10**(ls*(df['kic_steff']-5500)+li)]
        df = df.dropna(subset=['kic_smass']).reset_index()

    # Used for AS plots
    elif table=='cks+gaia2+h13':
        df1 = load_table('m17+gaia2+j17+iso')
        df1 = df1.groupby('id_kic',as_index=False).first()
        df2 = load_table('huber13')
        df = pd.merge(df1,df2)

    elif table=='cks+gaia2+s15':
        df1 = load_table('m17+gaia2+j17+iso')
        df1 = df1.groupby('id_kic',as_index=False).first()
        df2 = load_table('silva15')
        df = pd.merge(df1,df2)

    elif table=='m17+gaia2+j17+ext':
        df = load_table('m17+gaia2+j17')
        df['distance'] = np.array(1 / df.gaia2_sparallax * 1000) * u.pc
        df['ra'] = df['m17_ra']
        df['dec'] = df['m17_dec']
        df = cksgaia.extinction.add_extinction(df,'bayestar2017')
        df = df.drop('distance ra dec'.split(),axis=1)

    elif table=='m17+gaia2+j17+iso':
        df1 = load_table('m17+gaia2+j17+extinct')
        df2 = pd.read_csv(os.path.join(DATADIR, 'isochrones_gaia2.csv'))
        df = pd.merge(df1, df2, on='id_starname')

    elif table=='m17+gaia2+j17+iso+fur17':
        df1 = load_table('m17+gaia2+j17+iso')
        import pdb;pdb.set_trace()
        df2 = load_table('fur17')
        df = pd.merge(df1,df2 ,how='left')

    elif table == "cksgaia-planets":
        df2 = load_table('m17+gaia2+j17+iso+fur17')
        df = cksgaia.calc.update_planet_parameters(df2)

    elif table == "cksgaia-planets-filtered":
        df = load_table('cksgaia-planets')
        df = apply_filters(df)

    elif table == 'cksgaia-planets-weights':
        df = load_table('cksgaia-planets-filtered')
        df = cksgaia.completeness.weight_merge(df)

    else:
        assert False, "table {} not valid table name".format(table)
    return df

def read_silva15(fn):
    with open(fn,'r') as f:
        lines = f.readlines()
    header = lines[7]

    lines = lines[9:]
    _lines = []
    for line in lines:
        if line.count('&') > 0:
            _lines.append(line)

    lines = _lines

    _lines = []
    i=0
    df = []

    for line in lines:
        d = {}
        line =  line.split('&')
        d = OrderedDict()
        d['id_koi'] = line[0]
        d['id_kic'] = line[1]
        d['teff'] = line[2].split('$\\pm$')[0]
        d['teff_err1'] = line[2].split('$\\pm$')[1]

        d['fe'] = line[3].split('$\\pm$')[0]
        d['fe_err1'] = line[3].split('$\\pm$')[1]

        mass = re.sub(r"\$|\{|\^|\}|\_|\+"," ",line[4]).split()
        d['smass'] = mass[0]
        d['smass_err1'] = mass[1]
        d['smass_err2'] = mass[2]

        radius = re.sub(r"\$|\{|\^|\}|\_|\+"," ",line[5]).split()
        d['srad'] = radius[0]
        d['srad_err1'] = radius[1]
        d['srad_err2'] = radius[2]

        logg = re.sub(r"\$|\{|\^|\}|\_|\+"," ",line[7]).split()
        d['slogg'] = logg[0]
        d['slogg_err1'] = logg[1]
        d['slogg_err2'] = logg[2]

        age = re.sub(r"\$|\{|\^|\}|\_|\+"," ",line[9]).split()
        d['sage'] = age[0]
        d['sage_err1'] = age[1]
        d['sage_err2'] = age[2]

        df.append(d)
        i+=1

    df = pd.DataFrame(df).convert_objects(convert_numeric=True)
    df['teff_err2'] = -1.0 * df['teff_err1']
    df['slogage'] = np.log10(df['sage']) + 9 
    df['slogage_err1'] = np.log10(df.sage+df.sage_err1)+9 - df.slogage
    df['slogage_err2'] = np.log10(df.sage+df.sage_err2)+9 - df.slogage
    df['fe_err2'] = -1.0 * df['fe_err1']
    df = add_prefix(df,'s15_')
    return df
def read_huber13(fn, readme):
    df = astropy.io.ascii.read(fn, readme=readme)
    df = df.to_pandas()
    namemap = {
        'KIC':'id_kic',
        'KOI':'id_koi',
        'Mass':'h13_smass',
        'Rad':'h13_srad',
        'e_Mass':'h13_smass_err',
        'e_Rad':'h13_srad_err',
    }
    df = df.rename(columns=namemap)[namemap.values()]
    df = df.query('h13_srad > 0.5')
    df['id_kic'] = df.id_kic.astype(int)
    return df

def read_furlan17_table2(fn):
    df = pd.read_csv(fn,sep='\s+')
    namemap = {'KOI':'id_koi','KICID':'id_kic','Observatories':'ao_obs'}
    df = df.rename(columns=namemap)[namemap.values()]
    df['id_starname'] = ['K'+str(x).rjust(5, '0') for x in df.id_koi] 
    df = add_prefix(df,'fur17_')
    return df

def read_furlan17_table9(fn):
    names = """
    id_koi hst hst_err i i_err 692 692_err lp600 lp600_err jmag jmag_err 
    kmag kmag_err jkdwarf jkdwarf_err jkgiant jkgiant_err rcorr_avg 
    rcorr_avg_err
    """.split()

    df = pd.read_csv(fn,sep='\s+',skiprows=2,names=names)
    df['id_starname'] = ['K'+str(x).rjust(5, '0') for x in df.id_koi] 
    df = add_prefix(df,'fur17_')
    return df 

# Table manipulation

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

def order_columns(df, verbose=False, drop=True):
    columns = list(df.columns)
    coldefs = load_table('coldefs')
    cols = []
    for col in coldefs.column:
        if columns.count(col) == 1:
            cols.append(col)

    df = df[cols]
    if verbose and (len(cols) < len(columns)):
        print "table contains columns not defined in coldef"

    return df
