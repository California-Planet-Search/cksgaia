import pandas as pd

from matplotlib import pylab as plt
import numpy as np
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import cksgaia.io
import cksgaia
errorbar_kw = dict(markersize=6,color='b')

def comparison(key):
    if key=='smass-h13':
        df = cksgaia.io.load_table('cks+gaia2+h13')
        x1 = df.h13_smass
        x2 = df.giso_smass
        x3 = x2 / x1 

        x1err = np.vstack([df.h13_smass_err,df.h13_smass_err])
        x2err = np.vstack([-df.giso_smass_err2,df.giso_smass_err1])
        x3err = x2err / np.array(df.giso_smass)
        fig, axL = subplots_compare(
            x1,x2,x3, x1err=x1err, x2err=x2err, x3err=x3err, **errorbar_kw
        )
        _ylabel0 = '$M_\star$ [CKS] (Solar-masses)'
        _xlabel0 = '$M_\star$ [H13] (Solar-masses)'
        _ylabel1 = 'CKS / H13'
        _ylim1 = 0.8,1.2
        _xlim0 = 0.7,1.2
        _yt1 = [0.9,1.0,1.1]
        _xt0 = arange(0.8,2.0,0.2)

    if key=='srad-h13':
        df = cksgaia.io.load_table('cks+gaia2+h13')
        df = df.query('gaia2_gflux_ratio < 1.1 and gdir_srad_err1/gdir_srad < 0.1')

        x1 = df.h13_srad
        x2 = df.gdir_srad
        x3 = x2 / x1 
        x1err = np.vstack([df.h13_srad_err,df.h13_srad_err]) 
        x2err = np.vstack([-df.gdir_srad_err2,df.gdir_srad_err1])
        x3err = x2err / np.array(df.gdir_srad)
        fig, axL = subplots_compare(
            x1,x2,x3, x1err=x1err, x2err=x2err, x3err=x3err, **errorbar_kw
        )
        axL[0].set_xscale('log')
        axL[0].set_yscale('log')
        axL[1].set_xscale('log')    
        axL[1].set_yscale('linear')
        _ylabel0 = '$R_\star$ [CKS] (Solar-radii)'
        _xlabel0 = '$R_\star$ [H13] (Solar-radii)'
        _ylabel1 = 'CKS / H13'
        _ylim1 = 0.7,1.3
        _xlim0 = 0.4,12
        _yt1 = [0.8,0.9,1.0,1.1,1.2]
        _xt0 = [0.5,1,2,3,5,10,20]


    if key=='sage-s15':
        df = cksgaia.io.load_table('cks+gaia2+s15')
        x1 = df.s15_sage
        x2 = df.giso_sage
        x3 = x2 / x1 
        x1err = np.vstack([-df.s15_sage_err2,df.s15_sage_err1])
        x2err = np.vstack([-df.giso_sage_err2,df.giso_sage_err1])
        x3err = x2err / np.array(df.giso_sage)
        fig, axL = subplots_compare(
            x1,x2,x3, x1err=x1err, x2err=x2err, x3err=x3err, **errorbar_kw
        )
        axL[0].set_xscale('log')
        axL[0].set_yscale('log')
        axL[1].set_xscale('log')    
        axL[1].set_yscale('log')
        _ylabel0 = 'Age [CKS] (Gyr)'
        _xlabel0 = 'Age [S15] (Gyr)'
        _ylabel1 = 'CKS / S15'
        _ylim1 = 0.33,3
        _xlim0 = 1.5,15
        _yt1 = [0.5,1,2]
        _xt0 = [2,3,5,7,10,15]
        axL[0].minorticks_off()
        axL[1].minorticks_off()

class Comparison(object):
    """ 
    Code for generic comparisons (all parameters)
    """
    def __init__(self):
        pass

    def _provision_figure(self):
        self.fig = plt.figure(figsize=figsize)
        ax1 = plt.subplot2grid((4,1), (0,0), rowspan=3)
        ax2 = plt.subplot2grid((4,1), (3,0), rowspan=1, sharex=ax1)
        self.axL = [ax1,ax2]

    def _label_figure(self):
        plt.sca(self.axL[0])
        plt.setp(self.axL[0],ylabel=self.ylabel0,xlim=self.xlim0,ylim=self.xlim0)
        plt.xticks(self.xt0,self.xt0)
        plt.yticks(self.xt0,self.xt0)
        plt.sca(self.axL[1])
        plt.setp(self.axL[1],ylabel=self.ylabel1,xlabel=self.xlabel0,ylim=self.ylim1)
        plt.xticks(self.xt0,self.xt0)
        plt.yticks(self.yt1,self.yt1)

        one2one_kw = dict(linestyle='--',lw=1,color='g')
        plt.sca(self.axL[0])
        xt = plt.gca().get_xticklabels()
        plt.setp(xt,visible=False)

        one2one(**one2one_kw)
        plt.sca(self.axL[1])
        plt.axhline(1, **one2one_kw)

    def _subplots_compare(self, **errorbar_kw):
        subplots_compare(
            self.x1,self.x2,self.x3, self.fig, self.axL, 
            x1err=self.x1err, x2err=self.x2err, x3err=self.x3err, 
            **errorbar_kw
        )

class ComparisonRadius(Comparison):
    def __init__(self, key):
        if key=='srad-h13':
            df = cksgaia.io.load_table('cks+gaia2+h13')
            df = df.query('gaia2_gflux_ratio < 1.1 and gdir_srad_err1/gdir_srad < 0.1')
            
            x1 = df.h13_srad
            x2 = df.gdir_srad
            x3 = x2 / x1 
            x1err = np.vstack([df.h13_srad_err,df.h13_srad_err]) 
            x2err = np.vstack([-df.gdir_srad_err2,df.gdir_srad_err1])
            x3err = x2err / np.array(df.gdir_srad)
            xlabel0 = '$R_\star$ [H13] (Solar-radii)'
            ylabel0 = '$R_\star$ [CKS+Gaia] (Solar-radii)'

        if key=='srad-s15':
            df = cksgaia.io.load_table('cks+gaia2+s15')
            df = df.query('gaia2_gflux_ratio < 1.1 and gdir_srad_err1/gdir_srad < 0.1')
            x1 = df.s15_srad
            x2 = df.gdir_srad
            x3 = x2 / x1 
            x1err = np.vstack([-df.s15_srad_err2,df.s15_srad_err1]) 
            x2err = np.vstack([-df.gdir_srad_err2,df.gdir_srad_err1])
            x3err = x2err / np.array(df.gdir_srad)
            xlabel0 = '$R_\star$ [S15] (Solar-radii)'
            ylabel0 = '$R_\star$ [CKS+Gaia] (Solar-radii)'

        if key=='srad-j17':
            df = cksgaia.io.load_table('j17').groupby('id_kic',as_index=False).nth(0)
            df = df['id_kic iso_srad iso_srad_err1 iso_srad_err2'.split()]
            df = cksgaia.io.sub_prefix(df, 'iso', ignore=['id'])
            df = cksgaia.io.add_prefix(df, 'j17', ignore=['id'])
            df1 = df

            df2 = cksgaia.io.load_table('m17+gaia2+j17+iso+fur17')
            df2 = df2.groupby('id_kic',as_index=False).nth(0)
            df = pd.merge(df1,df2)
            df = df.query('gaia2_gflux_ratio < 1.1 and gdir_srad_err1/gdir_srad < 0.1 and ~(fur17_rcorr_avg > 1.05) and abs(gaia2_gmag - kic_kepmag) < 0.2')
#            df = df.query('gaia2_gflux_ratio < 1.1 and gdir_srad_err1/gdir_srad < 0.1 and ~(fur17_rcorr_avg > 1.05)')

            x1 = df.gdir_srad
            x2 = df.j17_srad
            x3 = x2 / x1 
            x1err = np.vstack([-df.gdir_srad_err2,df.gdir_srad_err1]) 
            x2err = np.vstack([-df.j17_srad_err2,df.j17_srad_err1])
            x3err = x2err / np.array(df.gdir_srad)
            xlabel0 = '$R_\star$ [CKS+Gaia] (Solar-radii)'
            ylabel0 = '$R_\star$ [CKS-II] (Solar-radii)'

        if key=='srad-gaia2':
            df = cksgaia.io.load_table('m17+gaia2+j17+iso+fur17')
            df = df.query('gaia2_gflux_ratio < 1.1 and gdir_srad_err1/gdir_srad < 0.1 and ~(fur17_rcorr_avg > 1.05)')
            df = df.groupby('id_kic',as_index=False).nth(0)

            x1 = df.gdir_srad
            x2 = df.gaia2_srad
            x3 = x2 / x1 
            x1err = np.vstack([-df.gdir_srad_err2,df.gdir_srad_err1]) 
            x2err = np.vstack([-df.gaia2_srad_err2,df.gaia2_srad_err1])
            x3err = x2err / np.array(df.gdir_srad)
            xlabel0 = '$R_\star$ [CKS+Gaia] (Solar-radii)'
            ylabel0 = '$R_\star$ [GaiaDR2] (Solar-radii)'

        self.df = df
        self.ylim1 = 0.7,1.3
        self.xlim0 = 0.4,12
        self.yt1 = [0.8,0.9,1.0,1.1,1.2]
        self.xt0 = [0.5,1,2,3,5,10,20]

        if (key.count('h13') + key.count('s15')) > 0:
            self.ylim1 = 0.85,1.15
            self.xlim0 = 0.4,12
            self.yt1 = [0.9,0.95,1.0,1.05,1.1]
            self.xt0 = [0.5,1,2,3,5,10,20] 

        self.x1 = x1
        self.x2 = x2
        self.x3 = x3 
        self.x1err = x1err
        self.x2err = x2err
        self.x3err = x3err
        self.ylabel0 = ylabel0
        self.xlabel0 = xlabel0
        self.ylabel1 = "Ratio $(y/x)$"

    def plot_comparison(self):
        self._provision_figure()
        self._subplots_compare()
        self.axL[0].set_xscale('log')
        self.axL[0].set_yscale('log')
        self.axL[1].set_xscale('log')    
        self.axL[1].set_yscale('linear')
        self._label_figure()
       
    def mean_string(self):
        _mean  = self.x3.mean() 
        if _mean > 1:
            s = r"{:.1f}\% larger".format(100*(_mean-1.0))
        else:
            s = r"{:.2f}\% smaller".format(abs(100*(_mean-1.0)))
        return s

    def std_string(self):
        std = self.x3.std()
        s = r"{:.1f}\%".format(100*std)
        return s

def subplots_compare(x1, x2, x3, fig, axL, x1err=None, x2err=None, x3err=None, **kwargs):
    mfc = 'RoyalBlue'
    fig.set_tight_layout(False)
    fig.subplots_adjust(hspace=0.001,left=0.20,top=0.95,right=0.90,bottom=0.13)
    plt.sca(axL[0])
    plt.plot(x1, x2, 'o',markersize=5,color='k')
    plt.errorbar(
        x1, x2, xerr=x1err, yerr=x2err, fmt='o',markersize=3.5,elinewidth=1,
        ecolor='DarkGray',color=mfc, capsize=0
    )

    plt.sca(axL[1])
    plt.plot(x1, x3, 'o',markersize=5,color='k')
    plt.errorbar(
        x1, x3, xerr=x1err, yerr=x3err, fmt='o',markersize=3.5,elinewidth=1,
        ecolor='DarkGray',color=mfc, capsize=0
    )
    return fig,axL

def add_anchored(*args,**kwargs):
    ax = plt.gca()
    at = AnchoredText(*args,**kwargs)
    ax.add_artist(at)

def one2one(**kwargs):
    xl = plt.xlim()
    plt.plot(xl,xl,**kwargs)

figsize = (4,4.5)
# figsize = (5,3.5)
