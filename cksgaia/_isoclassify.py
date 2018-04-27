import os 
import numpy as np
import pylab as pl
import pandas as pd
from matplotlib.pylab import * 
import time

import ebf

import grid.classify_grid 

from cksgaia.io import DATADIR
import cksgaia.iso
import cksgaia.extinction

def load_mist():
    model = ebf.read(os.path.join(DATADIR,'mesa.ebf'))
    return model

class Pipeline(cksgaia.iso.Pipeline):
    def run(self, dmodel=None):
        self.logg = -99
        self.logg_err = 99.0

        print "teff ", self.teff, self.teff_err
        print "logg ", self.logg, self.logg_err
        print "met ", self.met, self.met_err
        print "Kmag ", self.kmag, self.kmag_err
        print "par", self.parallax, self.parallax_err
        print "dist", 1/self.parallax, 1/self.parallax - 1/(self.parallax + self.parallax_err)

        model = load_mist()
        # prelims to manipulate some model variables (to be automated soon ...)
        model['rho']=np.log10(model['rho'])
        # next line turns off Dnu scaling relation corrections
        model['fdnu'][:]=1.
        model['avs']=np.zeros(len(model['teff']))
        model['dis']=np.zeros(len(model['teff']))


        # Instantiate model
        x = grid.classify_grid.obsdata()

        # add any combiantion of observables
        # Teff, logg, FeH + uncertainties
        x.addspec(
            [self.teff,self.logg,self.met],
            [self.teff_err,self.logg_err,self.met_err]
        )

        x.addplx(self.parallax, self.parallax_err)

        # Get extinction from bayestar model
        edf = pd.DataFrame([], columns=['ra', 'dec', 'parallax'])
        edf['ra'] = [self.ra]
        edf['dec'] = [self.dec]
        edf['gaia2_sparallax'] = [self.parallax]
        for i in range(10):
            try:
                edf = cksgaia.extinction.add_extinction(edf, 'bayestar2017')
                break
            except:
                time.sleep(2)
            if i == 10:
                print "WARNING: Extinction correction failed"

        self.kmag_ext = self.kmag + edf['ak'].values[0]
        self.kmag_ext_err = np.sqrt(self.kmag_err**2 + edf['ak_err'].values[0]**2)
        print "Kmag_ext ", self.kmag_ext, self.kmag_ext_err

        x.addjhk([-99,-99, self.kmag_ext],[0,0,self.kmag_ext_err])
        # Sloan photometry
        #x.addgriz([11.776,11.354,11.238,11.178],[0.02,0.02,0.02,0.02])
        paras = grid.classify_grid.classify(
            input=x, model=model, dustmodel=0, doplot=0, useav=0
        )
        
        #gcf('posteriors').savefig(self.pngfn)
        # fig = pl.figure('posteriors')
        # fig.savefig(self.pngfn)

        # fig = pl.figure('hrd')
        # fig.savefig(self.pngfn.replace('isochrones.png', 'hrd.png'))

        outdf = {}
        coldefs = {'iso_steff': 'teff',
                   'iso_slogg': 'logg',
                   'iso_smet': 'feh',
                   'iso_srad': 'rad',
                   'iso_smass': 'mass',
                   'iso_sage': 'age'}
        for outcol,incol in coldefs.items():
            outdf[outcol] = getattr(paras, incol)
            outdf[outcol+'_err1'] = getattr(paras, incol+'ep')
            outdf[outcol+'_err2'] = -getattr(paras, incol+'em')
        
        logage = np.log10(paras.age*1e9)
        logage_upper = np.log10(paras.age*1e9 + paras.ageep*1e9) - logage
        logage_lower = np.log10(paras.age*1e9 - paras.ageem*1e9) - logage
        outdf['iso_slogage'] = logage
        outdf['iso_slogage_err1'] = logage_upper
        outdf['iso_slogage_err2'] = logage_lower

        # par = 1 / paras.dis * 1e3
        # par_lower = 1 / (paras.dis + paras.disep) * 1e3 - par
        # par_upper = 1 / (paras.dis - paras.disem) * 1e3 - par
        # outdf['iso_sparallax'] = par
        # outdf['iso_sparallax_err1'] = par_upper
        # outdf['iso_sparallax_err2'] = par_lower
        outdf['iso_sparallax'] = self.parallax
        outdf['iso_sparallax_err1'] = self.parallax_err
        outdf['iso_sparallax_err2'] = -self.parallax_err
        outdf = pd.Series(outdf)
        outdf.to_csv(self.csvfn)
        
def _csv_reader(f):
    ifile = open(f, 'r')
    dirname = os.path.dirname(f)
    id_starname = dirname.split('/')[-1]
    logfile = open(f.replace('isochrones.csv', 'run-iso.log'), 'r')
    for line in logfile.readlines():
        if line.startswith("number of models after phot constraints"):
            num_mods = line.split()[-1]
            break
    row = {'id_starname': id_starname}
    for line in ifile.readlines():
        key, val = line.strip().split(',')
        key = key.replace('iso_', 'giso_')
        row[key] = float(val)
    row['giso_nmodels'] = int(num_mods)
    ifile.close()
    row = pd.DataFrame.from_dict(row, orient='index')
    row = row.T

    return row
