import cksgaia.io
#import seaborn as sns
from matplotlib.pylab import *
# sns.set_context('paper',font_scale=1.2)
# sns.set_style('ticks')

def fig_extinction():
    fig, axL = subplots(ncols=3,figsize=(7.5,2.5))
    sca(axL[0])
    fig_galactic_latitude()
    sca(axL[1])
    fig_distance()
    sca(axL[2])
    fig_compare()
    fig.set_tight_layout(True)

def fig_galactic_latitude():
    df = cksgaia.io.load_table('j17+m17+extinct',cache=1)
    semilogy()
    plot(df.m17_b,df.g15_a2massk,'.')
    xlabel('Galactic Latitude (deg)')
    ylabel('Ak (mag) [G15]')
    xlim(0,25)
    grid()
    
def fig_distance():
    df = cksgaia.io.load_table('j17+m17+extinct',cache=1)
    #df = df.query('6 < m17_b < 8 ')
    loglog()
    plot(df.sdistance,df.g15_a2massk,'.')
    xlabel('Isochrone Distance (pc)')
    ylabel('Ak (mag) [G15]')
    xlim(30,3000)
    grid()
    
def fig_compare():
    df = cksgaia.io.load_table('j17+m17+extinct',cache=1)
    plot(df.d03_a2massk,df.d03_a2massk-df.g15_a2massk,'.')
    xlim(0,0.15)
    ylim(-0.050,0.05)
    xlabel('Ak (mag) [D03]')
    ylabel('$\Delta$ Ak (mag) [D03 - G15]')
    xlim(0,0.2)
    ylim(-0.05,0.05)
    grid()


