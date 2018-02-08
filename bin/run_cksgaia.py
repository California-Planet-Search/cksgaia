#!/usr/bin/env python
from argparse import ArgumentParser
import os
from collections import OrderedDict

import glob
import cksgaia.io     # module for reading and writing datasets
import cksgaia.value  # module for computing scalar values for table
import cksgaia.table  # module for computing scalar values for table
import cksgaia.plot   # submodule for including plots


def main():
    psr = ArgumentParser()
    subpsr = psr.add_subparsers(title="subcommands", dest='subcommand')
    psr_parent = ArgumentParser(add_help=False)

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
        'isoclassify','isochrones'
    ]
    psr2.add_argument('mode',choices=modes)
    psr2.add_argument('baseoutdir')
    psr2.add_argument('outfile')
    psr2.set_defaults(func=create_iso_table)

    psr2 = subpsr.add_parser('create-val', parents=[psr_parent], )
    psr2.add_argument('name',type=str)
    psr2.set_defaults(func=create_val)

    psr2 = subpsr.add_parser('create-plot', parents=[psr_parent], )
    psr2.add_argument('name',type=str)
    psr2.set_defaults(func=create_plot)

    psr2 = subpsr.add_parser('create-table', parents=[psr_parent], )
    psr2.add_argument('name',type=str)
    psr2.set_defaults(func=create_table)

    psr2 = subpsr.add_parser('update-paper', parents=[psr_parent])
    psr2.set_defaults(func=update_paper)

    args = psr.parse_args()
    args.func(args)

def run_iso(args):
    import cksgaia.iso
    cksgaia.iso.run(args.driver, args.id_starname, args.outdir,debug=args.debug)

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
        print "mkdir -p {}; run_cksgaia.py run-iso {} {} {} &> {}/run-iso.log".format(outdir, args.driver, id_starname, outdir, outdir)
    
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


    




def create_table(args):
    w = Workflow()
    w.create_file('table', args.name ) 


def create_plot(args):
    w = Workflow()
    w.create_file('plot', args.name ) 


def create_val(args):
    w = Workflow()
    w.create_file('val',args.name) 


def update_paper(args):
    w = Workflow()
    w.update_paper()


class Workflow(object):
    def __init__(self):
        d = OrderedDict()

        # register different plots here
        # d['lamo-on-cks'] = cksmet.plotting.calibrate.lamo_on_cks
        # run_cksgaia create-plot lamo-on-cks # fig_lamo-on-cks.pdf, fig_lamo-on-cks.pdf fig_lamo-on-cks.pdf
        # run_cksgaia create-plot #

        self.plot_dict = d

        d = OrderedDict()
        # register different tables here
        # d['population-stub'] = lambda : cksmet.tables.population(stub=True)
        self.table_dict = d

        d = OrderedDict()
        #d['fit'] = cksmet.values.val_fit
        self.val_dict = d

        d = OrderedDict()
        d['table'] = self.table_dict
        d['plot'] = self.plot_dict
        d['val'] = self.val_dict
        self.all_dict = d

    def key2fn(self, key, kind):
        if kind=='plot':
            return 'fig_'+key+'.pdf'
        if kind=='table':
            return 'tab_'+key+'.tex'
        if kind=='val':
            return 'val_'+key+'.tex'
            
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
                gcf().savefig(fn)

            elif kind=='table':
                if name=='all':
                    lines = func()
                elif name==key:
                    lines = func()
                else:
                    continue
                    
                # Remove last \\
                fn = self.key2fn(key, 'table')
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

