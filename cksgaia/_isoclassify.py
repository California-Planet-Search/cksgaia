import os 
import numpy as np
import pylab as pl
import pandas as pd
from matplotlib.pylab import * 

import grid.classify_grid 

import cksgaia.io
import cksgaia.iso

import mwdust

class Pipeline(cksgaia.iso.Pipeline):
    def run(self, dmodel=None):

        print "teff ", self.teff, self.teff_err
        print "logg ", self.logg, self.logg_err
        print "met ", self.met, self.met_err
        print "Kmag ", self.kmag, self.kmag_err
        print "par", self.parallax, self.parallax_err
        print "dist", 1/self.parallax, 1/self.parallax - 1/(self.parallax + self.parallax_err)

        model = cksgaia.io.load_mist()
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
        edf['ra'] = self.ra
        edf['dec'] = self.dec

        x.addjhk([-99,-99, self.kmag],[0,0,self.kmag_err])
        # Sloan photometry
        #x.addgriz([11.776,11.354,11.238,11.178],[0.02,0.02,0.02,0.02])
        paras = grid.classify_grid.classify(
            input=x,model=model,dustmodel=0, doplot=0
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
    row = {'id_starname': id_starname}
    for line in ifile.readlines():
        key, val = line.strip().split(',')
        row[key] = float(val)
    ifile.close()
    row = pd.DataFrame.from_dict(row, orient='index')
    row = row.T
    return row
