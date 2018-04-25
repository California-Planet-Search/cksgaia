"""
Module for computing excintion to kepler targets
"""
import cksgaia.io
from astropy.coordinates import SkyCoord
from astropy import units as u
import os
import numpy as np
import pandas as pd
from dustmaps.bayestar import BayestarWebQuery

def add_extinction(df, mode):
    """
    add_extinction
    
    Args:
        df (pandas.DataFrame): must have following columns:
            - ra (decimal degrees)
            - dec (decimal degrees)
            - parallax (mas)
        mode (str): which dust model to use
            - bayestar2017 (Green et al. 2018)
            - bayestar2015 (Green et al. 2015)
    
    Returns:
        pandas.DataFrame: with following columns added
            - ak: extinction in K
            - ak_err: error on extinction in K including Rv and E(B-V)
    """
    dist = np.array(1 / df['gaia2_sparallax'] * 1000) * u.pc

    coords = SkyCoord(ra=np.array(df.ra)*u.degree, dec=np.array(df.dec)*u.degree, distance=dist, frame='icrs')
    rk_frac_err = 0.3 # Fractional uncertainty in R_K
    if mode=='bayestar2017':
        bayestar = BayestarWebQuery(version='bayestar2017')
        rk = 0.224 # A_K / E(B-V) 
        
    if mode=='bayestar2015':
        bayestar = BayestarWebQuery(version='bayestar2017')
        rk = 0.310 # A_K / E(B-V) 

    ebv = bayestar(coords, mode='percentile', pct=[16., 50., 84.])

    ak = rk * ebv

    ak_err1_map = ak[:,2] - ak[:,1] # formal A_K due to E(B-V)
    ak_err2_map = ak[:,0] - ak[:,1] # formal A_K due to E(B-V)
    ak_err_map = 0.5*(ak[:,2] - ak[:,0])

    ak_err_rk = ak[:,1] * rk_frac_err
    ak_err = ak[:,1] * np.sqrt( (ak_err_map / ak[:,1])**2 + rk_frac_err**2 )

    df['ak'] = ak[:,1]
    df['ak_err'] = ak_err
    df['ak_err_map'] = ak_err_map
    df['ebv'] = ebv[:,1]
    df['ebv_err'] = 0.5*(ebv[:,2] - ebv[:,0])
    df = pd.DataFrame(df)
    return df





