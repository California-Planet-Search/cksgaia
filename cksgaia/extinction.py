"""
Module for computing excintion to kepler targets
"""
import cksgaia.io
from astropy.coordinates import SkyCoord
from astropy import units as u
import os
import numpy as np
KEYS = [
    'd03_a2massk', # Drimmel (2003) 2MASS-K
    'd03_a2massj', # Drimmel (2003) 2MASS-J
    'd03_a2massh', # Drimmel (2003) 2MASS-H
    'd03_airac1', # Drimmel (2003) IRAC-1
    'd03_airac2', # Drimmel (2003) IRAC-2
    'g15_a2massk', # Green (2015) 
    'g15_a2massj', 
    'g15_a2massh',
    'g15_airac1', 
    'g15_airac2'
]

TABLES = [
    'j17+m17'
]

def compute(table, key, outdir='./', debug=False):
    print key
    smodel, sfilter = key.split('_') 

    if smodel=="g15":
        import mwdust.Green15 as model 
    elif smodel=="d03":
        import mwdust.Drimmel03 as model 
    else:
        assert False, "{} not supported".format(smodel)
        
    if sfilter=='a2massj':
        sfilter='2MASS J'
    elif sfilter=='a2massh':
        sfilter='2MASS H'
    elif sfilter=='a2massk':
        sfilter='2MASS Ks'
    elif sfilter=='airac1':
        sfilter='IRAC-1'
    elif sfilter=='airac2':
        sfilter='IRAC-2'

    if table=='j17+m17':
        df = cksgaia.io.load_table('j17+m17',cache=1)
        df = df.groupby('id_kic',as_index=False).first()
        df = df['id_kic m17_ra m17_dec iso_sparallax'.split()]
        df = df.rename(columns={'iso_sparallax':'sparallax'}) 

    print "computing {} band extinction using {} parallax and {} model".format(sfilter, table, smodel)

    mod = model(filter=sfilter)

    c = SkyCoord(ra=np.array(df.m17_ra)*u.degree, dec=df.m17_dec*u.degree, frame='icrs')
    df['m17_l'] = c.galactic.l.deg
    df['m17_b'] = c.galactic.b.deg
    df['sdistance'] = 1 / (df.sparallax * 1e-3)

    if debug:
        df = df.head()

    icount = 0 
    for i,row in df.iterrows():
        df.loc[i,key] = mod(row.m17_l,row.m17_b,row.sdistance)
        icount+=1
        if (icount % 100)==0:
            print "completed {}/{}".format(icount,len(df))
    
    fn = os.path.join(outdir, '{}-{}.csv'.format(table,key))
    print "saving to {}".format(fn)
    df.to_csv(fn)

    

