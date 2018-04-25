import seaborn as sns
import pandas as pd

from matplotlib.pylab import *
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import cksgaia.io
import cksgaia
errorbar_kw = dict(markersize=6,color='b')

def comparison(key):
    sns.set(style='ticks',rc={
        'ytick.major.size':3.0,'xtick.major.size':3.0,
        'xtick.direction': u'in','ytick.direction': u'in',}
    )

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
        df = df.query('gaia2_gflux_ratio < 1.1 and giso_srad_err1/giso_srad < 0.1')

        x1 = df.h13_srad
        x2 = df.giso_srad
        x3 = x2 / x1 
        x1err = np.vstack([df.h13_srad_err,df.h13_srad_err]) 
        x2err = np.vstack([-df.giso_srad_err2,df.giso_srad_err1])
        x3err = x2err / np.array(df.giso_srad)
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

    if key=='srad-j17':
        df = cksgaia.io.load_table('j17').groupby('id_kic',as_index=False).nth(0)
        df = df['id_kic iso_srad iso_srad_err1 iso_srad_err2'.split()]
        df = cksgaia.io.sub_prefix(df, 'iso', ignore=['id'])
        df = cksgaia.io.add_prefix(df, 'j17', ignore=['id'])
        df1 = df

        df2 = cksgaia.io.load_table('j17+m17+gaia2+iso+fur17').groupby('id_kic',as_index=False).nth(0)
        df = pd.merge(df1,df2)
        df = df.query('gaia2_gflux_ratio < 1.1 and giso_srad_err1/giso_srad < 0.1 and ~(fur17_rcorr_avg > 1.05)')

        x1 = df.giso_srad
        x2 = df.j17_srad
        x3 = x2 / x1 
        
        x1err = np.vstack([-df.giso_srad_err2,df.giso_srad_err1]) 
        x2err = np.vstack([-df.j17_srad_err2,df.j17_srad_err1])
        x3err = x2err / np.array(df.giso_srad)
        fig, axL = subplots_compare(
            x1,x2,x3, x1err=x1err, x2err=x2err, x3err=x3err, **errorbar_kw
        )
        axL[0].set_xscale('log')
        axL[0].set_yscale('log')
        axL[1].set_xscale('log')    
        axL[1].set_yscale('linear')
        _ylabel0 = '$R_\star$ [CKS-II] (Solar-radii)'
        _xlabel0 = '$R_\star$ [CKS+Gaia] (Solar-radii)'
        _ylabel1 = 'CKS / S15'
        _ylim1 = 0.7,1.3
        _xlim0 = 0.4,12
        _yt1 = [0.8,0.9,1.0,1.1,1.2]
        _xt0 = [0.5,1,2,3,5,10,20]
        import pdb;pdb.set_trace()


    if key=='srad-s15':
        df = cksgaia.io.load_table('cks+gaia2+s15')
        df = df.query('gaia2_gflux_ratio < 1.1 and giso_srad_err1/giso_srad < 0.1')
        x1 = df.s15_srad
        x2 = df.giso_srad
        x3 = x2 / x1 
        x1err = np.vstack([-df.s15_srad_err2,df.s15_srad_err1]) 
        x2err = np.vstack([-df.giso_srad_err2,df.giso_srad_err1])
        x3err = x2err / np.array(df.giso_srad)
        fig, axL = subplots_compare(
            x1,x2,x3, x1err=x1err, x2err=x2err, x3err=x3err, **errorbar_kw
        )
        axL[0].set_xscale('log')
        axL[0].set_yscale('log')
        axL[1].set_xscale('log')    
        axL[1].set_yscale('linear')
        _ylabel0 = '$R_\star$ [CKS] (Solar-radii)'
        _xlabel0 = '$R_\star$ [S15] (Solar-radii)'
        _ylabel1 = 'CKS / S15'
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

    sca(axL[0])
    setp(axL[0],ylabel=_ylabel0,xlim=_xlim0,ylim=_xlim0)
    xticks(_xt0,_xt0)
    yticks(_xt0,_xt0)
    sca(axL[1])
    setp(axL[1],ylabel=_ylabel1,xlabel=_xlabel0,ylim=_ylim1)
    xticks(_xt0,_xt0)
    yticks(_yt1,_yt1)
   
    one2one_kw = dict(linestyle='--',lw=1,color='g')
    sca(axL[0])
    xt = gca().get_xticklabels()
    setp(xt,visible=False)
    
    one2one(**one2one_kw)
    sca(axL[1])
    axhline(1, **one2one_kw)


def provision_figure():
    fig = figure(figsize=figsize)
    ax1 = subplot2grid((4,1), (0,0), rowspan=3)
    ax2 = subplot2grid((4,1), (3,0), rowspan=1, sharex=ax1)
    axL = [ax1,ax2]
    return fig,axL 

def subplots_compare(x1, x2, x3, xerr3=None, fig0=None, axL0=None, **kwargs):
    if fig0 is None:
        fig, axL = provision_figure()
    else:
        fig = fig0
        axL = axL0

    fig.set_tight_layout(False)
    fig.subplots_adjust(hspace=0.4,left=0.17,top=0.95,right=0.90)
    sca(axL[0])
    grid()
    errorbar(x1,x2,**kwargs)
    xl = xlim()
    yl = ylim()
    x12 = np.hstack([x1,x2])
    x = np.linspace(min(x12)/100,max(x12)*100)
    plot(x,x,zorder=1)
    xlim(*xl)
    ylim(*xl)

    sca(axL[1])
    grid() 
    errorbar(x1,x3,**kwargs)
    print "mean(x3) {}".format(np.mean(x3))
    print "std(x3) {}".format(np.std(x3))
    return fig,axL

def add_anchored(*args,**kwargs):
    ax = gca()
    at = AnchoredText(*args,**kwargs)
    ax.add_artist(at)

def provision_figure():
    fig = figure(figsize=figsize)
    ax1 = subplot2grid((4,1), (0,0), rowspan=3)
    ax2 = subplot2grid((4,1), (3,0), rowspan=1, sharex=ax1)
    axL = [ax1,ax2]
    return fig, axL 

def one2one(**kwargs):
    xl = xlim()
    plot(xl,xl,**kwargs)

figsize = (4,4.5)


def subplots_compare(x1, x2, x3, x1err=None, x2err=None, x3err=None, fig0=None, axL0=None, **kwargs):
    if fig0 is None:
        fig, axL = provision_figure()
    else:
        fig = fig0
        axL = axL0

    fig.set_tight_layout(False)

    #fig.subplots_adjust(hspace=0.4,left=0.17,top=0.95,right=0.90)
    fig.subplots_adjust(hspace=0.001,left=0.17,top=0.95,right=0.90,bottom=0.12)
    sca(axL[0])
    plot(x1, x2, 'o',markersize=5,color='k')
    errorbar(
        x1, x2, xerr=x1err, yerr=x2err, fmt='o',markersize=3.5,elinewidth=1,
        ecolor='DarkGray'
    )

    sca(axL[1])
    plot(x1, x3, 'o',markersize=5,color='k')
    errorbar(
        x1, x3, xerr=x1err, yerr=x3err, fmt='o',markersize=3.5,elinewidth=1,
        ecolor='DarkGray'
    )
    return fig,axL

def radius():
    df = cksphys.io.load_table('cks+iso+huber')
    x1 = df.huber_srad
    x2 = df.iso_srad
    x3 = x2 / x1 

    xerr = np.vstack([df.huber_srad_err,df.huber_srad_err])
    yerr = np.vstack([df.iso_srad_err1,-df.iso_srad_err2])

    fig, axL = provision_figure()
    axL[0].set_xscale('log')
    axL[0].set_yscale('log')
    axL[1].set_xscale('log')    
    axL[1].set_yscale('linear')


    sca(axL[0])
    grid()
    errorbar(x1,x2)
    xl = xlim()
    yl = ylim()
    x12 = np.hstack([x1,x2])
    x = np.linspace(min(x12)/100,max(x12)*100)
    plot(x,x,zorder=1)
    xlim(*xl)
    ylim(*xl)

    sca(axL[1])
    grid() 
    errorbar(x1,x3,fmt='.')
    print "mean(x3) {}".format(np.mean(x3))
    print "std(x3) {}".format(np.std(x3))


