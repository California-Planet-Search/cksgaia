#!/usr/bin/env python
from argparse import ArgumentParser
import os
from collections import OrderedDict
import glob

import pylab as pl
import pandas as pd

import cksgaia.io # module for reading and writing datasets
import cksgaia.value # module for computing scalar values for table
import cksgaia.table # module for computing scalar values for table
import cksgaia.plot.sample 
import cksgaia.plot.contour
import cksgaia.plot.occur
import cksgaia.plot.sim
import cksgaia.errors
import cksgaia.calc
import cksgaia.config
import cksgaia.sim.simulations
from cksgaia.plot.compare import ComparisonRadius as CR

def main():
    psr = ArgumentParser()
    subpsr = psr.add_subparsers(title="subcommands", dest='subcommand')
    psr_parent = ArgumentParser(add_help=False)
    psr_parent.add_argument('-d', dest='outputdir', type=str, help="put files in this directory")

    psr2 = subpsr.add_parser('create-iso-batch', parents=[psr_parent],)
    psr2.set_defaults(func=create_iso_batch)

    psr2 = subpsr.add_parser('create-xmatch-table', parents=[psr_parent])
    psr2.set_defaults(func=create_xmatch_table)

    psr2 = subpsr.add_parser('create-extinction-jobs', parents=[psr_parent])
    psr2.set_defaults(func=create_extinction_jobs)

    psr2 = subpsr.add_parser('compute-extinction', parents=[psr_parent])
    psr2.add_argument('table', help='name of star')
    psr2.add_argument('key')
    psr2.set_defaults(func=compute_extinction)

    psr2 = subpsr.add_parser('simulate-surveys', parents=[psr_parent])
    psr2.set_defaults(func=sim_surveys)

    psr_merge = subpsr.add_parser(
        'create-merged-table', parents=[psr_parent],
        description="Generate merged table with all of the columns."
    )
    psr_merge.set_defaults(func=create_merged_table)

    psr_stats = subpsr.add_parser(
        'tex-stats', parents=[psr_parent],
        description="Generate merged table with all of the columns."
    )
    psr_stats.set_defaults(func=tex_stats)

    psr2 = subpsr.add_parser('create-val', parents=[psr_parent], )
    psr2.add_argument('name',type=str)
    psr2.set_defaults(func=create_val)

    psr2 = subpsr.add_parser('create-plot', parents=[psr_parent], )
    psr2.add_argument('name',type=str)
    psr2.set_defaults(func=create_plot)

    psr2 = subpsr.add_parser('create-val', parents=[psr_parent], )
    psr2.add_argument('name',type=str)
    psr2.set_defaults(func=create_val)

    psr2 = subpsr.add_parser('create-table', parents=[psr_parent], )
    psr2.add_argument('name',type=str)
    psr2.set_defaults(func=create_table)

    psr2 = subpsr.add_parser('create-csv', parents=[psr_parent], )
    psr2.add_argument('name',type=str)
    psr2.set_defaults(func=create_csv)


    psr2 = subpsr.add_parser('update-paper', parents=[psr_parent])
    psr2.set_defaults(func=update_paper)

    args = psr.parse_args()
    args.func(args)


def create_xmatch_table(args):
    cksgaia.xmatch.create_xmatch_table()

def create_iso_batch(args):
    """Create Isoclassify Batch Jobs

    Creates input parameters for two runs
    
       1. The direct method with the following constraints
          - teff, logg, fe, parallax, kmag
       2. The grid method with the following constraints
          - teff, logg, met, kmag [no parallax]

    We default to the direct method. But if the parallax method from
    the grid based method is significantly different than the gaia
    parallax, there is additional flux in the aperture which indicates
    dilution.

    """
    df = cksgaia.io.load_table('m17+gaia2+j17').groupby('id_kic').nth(0)
    df = df.sort_values(by='id_starname')

    # Direct method
    df = df.rename(
        columns={
            'id_name':'id_starname',
            'parallax_error':'parallax_err',
            'ks_m':'kmag',
            'ks_msigcom':'kmag_err',
            'cks_steff':'teff',
            'cks_steff_err1':'teff_err',
            'cks_slogg':'logg',
            'cks_slogg_err1':'logg_err',
            'cks_smet':'feh',
            'cks_smet_err1':'feh_err',
            'm17_kmag':'kmag',
            'm17_kmag_err':'kmag_err',
            'gaia2_sparallax':'parallax',
            'gaia2_sparallax_err':'parallax_err',
            'gaia2_ra':'ra',
            'gaia2_dec':'dec',
        }
    )
    df['band'] = 'kmag'
    df['teff_err'] = 60
    df['parallax'] /= 1e3
    df['parallax_err'] /= 1e3
    df['feh_err'] = 0.04
    df.id_starname = df.id_starname.str.replace(' ','_')
    fn = 'data/isoclassify-direct.csv'
    df.to_csv(fn)
    print "created {}".format(fn)

    # Grid method. Don't set parallax so we can compare later
    df['parallax'] = -99
    df['parallax_err'] = 0
    fn = 'data/isoclassify-grid.csv'
    df.to_csv(fn)
    print "created {}".format(fn)


def sim_surveys(args):
    cksgaia.sim.simulations.run(args)

def create_extinction_jobs(args):
    for table in cksgaia.extinction.TABLES:
        for key in cksgaia.extinction.KEYS:
            print "run_cksgaia.py compute-extinction {} {}".format(table,key)

def compute_extinction(args):
    outdir = os.path.join(cksgaia.io.DATADIR,'extinction/')
    cmd = 'mkdir -p {}'.format(outdir)
    print cmd
    os.system(cmd)
    cksgaia.extinction.compute(args.table,args.key,outdir=outdir)

def create_iso_table(args):
    """
    Read in isochrones csvfiles 
    Args:
        outdir (str): where to look for isochrones.csv files
    """
    fL = glob.glob("{}/*/*.csv".format(args.baseoutdir))
    df = []

    import cksgaia._isoclassify
        
    if args.mode=='isoclassify':
        _csv_reader = cksgaia._isoclassify._csv_reader
    elif args.mode=='isochrones':
        _csv_reader = cksgaia._isochrones._csv_reader
    else:
        assert False, "invalid mode"

    for i, f in enumerate(fL):
        if i%100==0:
            print i
            
        try:
            df.append(_csv_reader(f))
        except ValueError:
            print "{} failed".format(f)

    df = pd.concat(df)
    df.to_csv(args.outfile, index=False)
    print "created {}".format(args.outfile)

def create_merged_table(args):
    df = cksgaia.io.load_table('j17+m17+gaia2+iso', verbose=True, cache=0)
    csvfn = os.path.join(cksgaia.io.DATADIR, 'cks_iso_gaia2_merged.csv')
    df.to_csv(csvfn)

def tex_stats(args):
    cksgaia.calc.table_statistics()

def create_table(args):
    w = Workflow(outputdir=args.outputdir)
    w.create_file('table', args.name)

def create_csv(args):
    w = Workflow(outputdir=args.outputdir)
    w.create_file('csv', args.name)

def create_plot(args):
    w = Workflow(outputdir=args.outputdir)
    w.create_file('plot', args.name)

def create_val(args):
    w = Workflow(outputdir=args.outputdir)
    w.create_file('val', args.name)

def update_paper(args):
    w = Workflow(outputdir=args.outputdir)
    w.update_paper() 

class Workflow(object):
    def __init__(self, outputdir='./'):
        self.outputdir = outputdir

        d = OrderedDict()

        # register different plots here
        d['insol-contour-masscuts'] = cksgaia.plot.contour.contour_masscuts
        d['period-contour-masscuts'] = cksgaia.plot.contour.period_contour_masscuts
        # d['extinction'] = cksgaia.plot.extinction.fig_extinction
        # d['sample'] = cksgaia.plot.sample.hrplot
        # d['filters'] = cksgaia.plot.sample.filter_plot
        # d['mag-hist'] = cksgaia.plot.sample.magcuts
        # d['depth-hist'] = cksgaia.plot.sample.depth_hist
        d['srad-h13'] = lambda: CR('srad-h13').plot_comparison()
        d['srad-s15'] = lambda: CR('srad-s15').plot_comparison()
        d['srad-j17'] = lambda: CR('srad-j17').plot_comparison()
        d['srad-gaia2'] = lambda: CR('srad-gaia2').plot_comparison()
        # d['smass-h13'] = lambda : CR('smass-h13').plot_comparison()
        # d['sage-s15'] = lambda : CR('sage-s15').plot_comparison()
        # d['srad-hist'] = cksgaia.plot.sample.srad_hist
        d['srad-err-hist'] = cksgaia.plot.sample.srad_err_hist
        d['prad-err-hist'] = cksgaia.plot.sample.prad_err_hist
        # d['parallax-err-hist'] = cksgaia.plot.sample.parallax_err_hist
        # d['insol-hist'] = cksgaia.plot.occur.insol_hist
        # d['radius-hist-fit'] = cksgaia.plot.occur.money_plot_fit
        d['radius-hist-plain'] = cksgaia.plot.occur.money_plot_plain
        # d['period-contour-q16'] = cksgaia.plot.contour.period_contour_q16
        d['period-contour-cks'] = cksgaia.plot.contour.period_contour_cks
        # d['insol-contour-anno'] = cksgaia.plot.contour.insol_contour_anno
        d['insol-contour-data'] = cksgaia.plot.contour.insol_contour_data
        # d['srad-contour'] = cksgaia.plot.contour.srad_contour
        # d['smass-cuts'] = cksgaia.plot.occur.mass_cuts
        # d['desert-edge'] = cksgaia.plot.occur.desert_edge
        d['desert-edge-cum'] = cksgaia.plot.occur.desert_edge_cum
        d['mean-values'] = cksgaia.plot.occur.mean_values
        d['width-sim-plot'] = cksgaia.plot.sim.wid_sim_plot
        d['per-prad'] = cksgaia.plot.sample.fig_per_prad
        d['insol-prad'] = cksgaia.plot.sample.fig_insol_prad
        d['mean-met'] = cksgaia.plot.occur.mean_met

        self.plot_dict = d

        # register different tables here
        d = OrderedDict()
        d['weight-tex-stub'] = lambda : cksgaia.table.weight_table(lines=10)
        d['weight-tex-full'] = lambda: cksgaia.table.weight_table(lines='all')
        d['star-stub'] = lambda: cksgaia.table.star()[:12]
        d['planet-stub'] = lambda: cksgaia.table.planet()[:15]
        #d['histbins'] = lambda: cksgaia.table.bins_table()
        #d['filters'] = lambda: cksgaia.table.filters_table()
        self.table_dict = d

        d = OrderedDict()
        # register machine-readable tables here
        d['weight-machine'] = cksgaia.table.weight_table_machine
        d['star-machine'] = cksgaia.table.star_machine
        d['planet-machine'] = cksgaia.table.planet_machine
        self.csv_dict = d

        d = OrderedDict()
        d['stat'] = cksgaia.value.val_stat
        self.val_dict = d

        d = OrderedDict()
        d['table'] = self.table_dict
        d['csv'] = self.csv_dict
        d['plot'] = self.plot_dict
        d['val'] = self.val_dict
        self.all_dict = d

    def key2fn(self, key, kind):
        if kind=='plot':
            return os.path.join(self.outputdir, 'fig_'+key+'.pdf')
        if kind=='table':
            return os.path.join(self.outputdir, 'tab_'+key+'.tex')
        if kind=='csv':
            return os.path.join(self.outputdir, 'tab_'+key+'.csv')
        if kind=='val':
            return os.path.join(self.outputdir, 'val_'+key+'.tex')
            
    def create_file(self, kind, name):
        i = 0
        for key, func in self.all_dict[kind].iteritems():
            if kind=='plot':
                if name=='all':
                    func()
                elif name==key:
                    func()
                else:
                    continue
                    
                fn = self.key2fn(key, 'plot')
                pl.gcf().savefig(fn)
                pl.clf()
                fig = pl.figure()

            elif kind=='table' or kind == 'csv':
                if name=='all':
                    lines = func()
                elif name==key:
                    lines = func()
                else:
                    continue
                    
                # Remove last \\
                fn = self.key2fn(key, kind)
                with open(fn,'w') as f:
                    f.writelines("\n".join(lines))

            elif kind=='val':
                fn = self.key2fn(key, 'val')
                if name=='all':
                    lines = func()
                elif name==key:
                    lines = func()
                else:
                    continue

                lines1 = [
                    "\\newcommand{\%s}[1]{%%" % key,
                    "\IfEqCase{#1}{",
                ]

                lines2 = [
                    "}[XX]",
                    "}"
                ]
                lines = lines1 + lines + lines2

                with open(fn,'w') as f:
                    f.writelines("%\n".join(lines))

            i+=1

        if i==0:
            assert False, name + " not a valid key"

    def update_paper(self):
        for kind, d in self.all_dict.iteritems():
            for key, val in d.iteritems():
                fn = self.key2fn(key, kind)
                cmd = 'cp {} paper/'.format(fn)
                print cmd
                os.system(cmd)

if __name__=="__main__":
    main()

