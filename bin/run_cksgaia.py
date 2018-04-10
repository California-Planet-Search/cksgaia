#!/usr/bin/env python
from argparse import ArgumentParser
import os
from collections import OrderedDict

import pylab as pl
import numpy as np

import pandas as pd
import glob

import cksgaia.io     # module for reading and writing datasets
import cksgaia.value  # module for computing scalar values for table
import cksgaia.table  # module for computing scalar values for table
import cksgaia.plot   # submodule for including plots
import cksgaia.errors
import cksgaia.calc
import cksgaia.plot.extinction
import cksgaia.extinction

def main():
    psr = ArgumentParser()
    subpsr = psr.add_subparsers(title="subcommands", dest='subcommand')
    psr_parent = ArgumentParser(add_help=False)
    psr_parent.add_argument('-d', dest='outputdir', type=str, help="put files in this directory")

    psr2 = subpsr.add_parser(
        'create-iso-jobs', parents=[psr_parent], 
        description="Create isochrone processing batch jobs for all KOIs"
    )
    psr2.add_argument('driver')
    psr2.add_argument('sample')
    psr2.add_argument('baseoutdir',help='absolute path to output directory')
    psr2.set_defaults(func=create_iso_jobs)

    psr2 = subpsr.add_parser(
        'run-iso', parents=[psr_parent], 
        description="Run isochrones"
    )

    drivers = [
        'isoclassify','isochrones','isocla+isochr-dsep','isocla+isochr-mist',
        'isocla+isochr-dsep-jk','isocla+isochr-mist-jk'
    ]
    psr2.add_argument('driver', help='isochrone interp code',choices=drivers)
    psr2.add_argument('id_starname', help='name of star')
    psr2.add_argument('outdir')
    psr2.add_argument('--debug',action='store_true')
    psr2.set_defaults(func=run_iso)

    psr2 = subpsr.add_parser(
        'create-iso-table', parents=[psr_parent], 
        description="Scrape the isochrones.csv files to make isochrones table"
    )

    modes = [
        'isoclassify', 'isochrones'
    ]
    psr2.add_argument('mode',choices=modes)
    psr2.add_argument('baseoutdir')
    psr2.add_argument('outfile')
    psr2.set_defaults(func=create_iso_table)

    psr2 = subpsr.add_parser('create-xmatch-table', parents=[psr_parent])
    psr2.set_defaults(func=create_xmatch_table)

    psr2 = subpsr.add_parser('create-extinction-jobs', parents=[psr_parent])
    psr2.set_defaults(func=create_extinction_jobs)

    psr2 = subpsr.add_parser('compute-extinction', parents=[psr_parent])
    psr2.add_argument('table', help='name of star')
    psr2.add_argument('key')
    psr2.set_defaults(func=compute_extinction)

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

    psr2 = subpsr.add_parser('update-paper', parents=[psr_parent])
    psr2.set_defaults(func=update_paper)

    args = psr.parse_args()
    args.func(args)

def run_iso(args):
    import cksgaia.iso
    cksgaia.iso.run(args.driver, args.id_starname, args.outdir, debug=args.debug)

def create_xmatch_table(args):
    cksgaia.io.create_xmatch_table()

def create_iso_jobs(args):
    if args.sample=='cks':
        df = cksgaia.io.load_table('j17+m17')
    else:
        import cksspec.io
        df = cksspec.io.load_table(args.sample)
        df['id_starname'] = df['name']

    for i, row in df.iterrows():
        id_starname = row.id_starname
        outdir = "{}/{}".format(args.baseoutdir, id_starname)
        print "mkdir -p {}; run_cksgaia.py run-iso {} {} {} &> {}/run-iso.log".format(outdir,
                                                                                      args.driver,
                                                                                      id_starname,
                                                                                      outdir, outdir)

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
        except:
            print "{} failed".format(f)

    df = pd.concat(df)
    df.to_csv(args.outfile, index=False)
    print "created {}".format(args.outfile)


def create_merged_table(args):
    df = cksgaia.io.load_table('cks+nea+iso', verbose=True, cache=0)

    csvfn = os.path.join(cksgaia.io.DATADIR, 'cks_fakegaia_merged.csv')
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
        d['extinction'] = cksgaia.plot.extinction.fig_extinction
        d['sample'] = cksgaia.plot.sample.hrplot
        d['filters'] = cksgaia.plot.sample.filter_plot
        d['mag-hist'] = cksgaia.plot.sample.magcuts
        d['depth-hist'] = cksgaia.plot.sample.depth_hist
        d['srad-hist'] = cksgaia.plot.sample.srad_hist
        d['srad-err-hist'] = cksgaia.plot.sample.srad_err_hist
        d['prad-err-hist'] = cksgaia.plot.sample.prad_err_hist
        d['insol-hist'] = cksgaia.plot.occur.insol_hist
        d['radius-hist-fit'] = cksgaia.plot.occur.money_plot_fit
        d['radius-hist-plain'] = cksgaia.plot.occur.money_plot_plain
        d['period-contour-q16'] = cksgaia.plot.contour.period_contour_q16
        d['period-contour-cks'] = cksgaia.plot.contour.period_contour_cks
        d['insol-contour-anno'] = cksgaia.plot.contour.insol_contour_anno
        d['insol-contour-data'] = cksgaia.plot.contour.insol_contour_data
        d['insol-contour-masscuts'] = cksgaia.plot.contour.contour_masscuts
        d['srad-contour'] = cksgaia.plot.contour.srad_contour
        d['smass-cuts'] = cksgaia.plot.occur.mass_cuts
        d['desert-edge'] = cksgaia.plot.occur.desert_edge
        d['desert-edge-cum'] = cksgaia.plot.occur.desert_edge_cum


        self.plot_dict = d

        d = OrderedDict()
        # register different tables here
        d['weight-tex-stub'] = lambda : cksgaia.table.weight_table(lines=10)
        d['weight-tex-full'] = lambda: cksgaia.table.weight_table(lines='all')
        d['histbins'] = lambda: cksgaia.table.bins_table()
        d['filters'] = lambda: cksgaia.table.filters_table()
        self.table_dict = d

        d = OrderedDict()
        # register machine-readable tables here
        d['weight-machine'] = lambda: cksgaia.table.weight_table_machine()
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
                    "\IfEqCase{#1}{%",
                ]

                lines2 = [
                    "}[\PackageError{tree}{Undefined option to tree: #1}{}]%",
                    "}%"
                ]
                lines = lines1 + lines + lines2

                with open(fn,'w') as f:
                    f.writelines("\n".join(lines))

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

