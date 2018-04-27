import glob
import os 

import numpy as np    
import pandas as pd
import cksgaia.io
import cksgaia.config
np.random.seed(0) 

class Pipeline(object):
    def __init__(self, id_starname, outdir):
        self.id_starname = id_starname
        self.outdir = outdir

        df = cksgaia.io.load_table('m17+gaia2+j17', cache=0)
        g = df.groupby('id_starname',as_index=True)
        df = g.nth(0)
        star = df.ix[id_starname]
        self.teff = star.cks_steff
        self.teff_err = star.cks_steff_err1
        self.logg = star.cks_slogg
        self.logg_err = star.cks_slogg_err1
        self.met = star.cks_smet
        self.met_err = star.cks_smet_err1
        self.kmag = star.kic_kmag
        self.kmag_err = 0.02
        self.jmag = star.kic_jmag
        self.jmag_err = 0.02
        self.parallax = star['gaia2_sparallax'] / 1e3
        self.parallax_err = star['gaia2_sparallax_err'] / 1e3
        self.ra = star.m17_ra
        self.dec = star.m17_dec

        self.pngfn = os.path.join(outdir,'isochrones.png')
        self.csvfn = os.path.join(outdir,'isochrones.csv')

def run(driver, id_starname, outdir,debug=False):
    if driver.count('isochrones')==1:
        import cksgaia._isochrones
        pipe = cksgaia._isochrones.Pipeline(id_starname, outdir) 
        pipe.run(100,100)

    elif driver.count('isoclassify')==1:
        import cksgaia._isoclassify
        pipe = cksgaia._isoclassify.Pipeline(id_starname, outdir) 
        pipe.run()

    elif (driver.count('isocla+isochr-dsep')==1 or 
          driver.count('isocla+isochr-mist')==1 ):
        import cksgaia._isoclassify

        bands = 'k'
        if driver.count('jk')==1:
            bands = 'jk'

        nwalk = 1000
        nburn = 300 
        niter = 300 

        if debug:
            nwalk = 100
            nburn = 100
            niter = 100 

        if driver.count('mist')>0:
            model='mist'
        elif driver.count('dsep')>0:
            model='dsep'
        else:
            assert False, "Problem"


        print "First run isoclassify to get the approximate solution "
        pipe = cksgaia._isoclassify.Pipeline(id_starname, outdir) 
        pipe.run()

        #p0 += m0, age0, feh0, d0, AV0
        # m0 default is solar with a std of 0.3dex
        # age 9 (1gyr) spread of 0.5 dex 
        # feh0 (-2.5--0.5dex) spread of 0.5 dex 
        # d0 (some r^2 prior from 1--30,000 pc) 
        # AV0 (1-0) 
        row = pd.read_csv(pipe.csvfn,squeeze=True,header=None, index_col=0)
        dist = 1e3 / row.gaia2_sparallax # mas -> pc
        dist_err = dist * row.gaia2_sparallax_err1 / row.gaia2_sparallax
        slogage = np.log10(row.iso_sage * 1e9) # Gyr -> logage
        slogage = np.log10(row.iso_sage * 1e9) # Gyr -> logage

        p0 = np.zeros((nwalk,5))
        
        mass = row.iso_smass + np.random.randn(nwalk) * row.iso_smass_err1*0.1
        slogage = slogage + np.random.randn(nwalk) * 0.1
        smet = row.iso_smet + np.random.randn(nwalk) * 0.01
        dist = dist + np.random.randn(nwalk) * dist_err
        av = np.random.uniform(0,1,nwalk)

        p0[:,0] = mass
        p0[:,1] = slogage
        p0[:,2] = smet
        p0[:,3] = dist
        p0[:,4] = av

        pipe = cksgaia._isochrones.Pipeline(id_starname, outdir) 
        pipe.run(nburn, niter, nwalk, p0, model, bands)

    else:
        assert False, "{} not a valid driver".format(driver)





