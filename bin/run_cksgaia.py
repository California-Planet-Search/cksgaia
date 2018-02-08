#!/usr/bin/env python
from argparse import ArgumentParser
import os
from collections import OrderedDict

import cksgaia.io     # module for reading and writing datasets
import cksgaia.value  # module for computing scalar values for table
import cksgaia.table  # module for computing scalar values for table
import cksgaia.plot   # submodule for including plots


def main():
    psr = ArgumentParser()
    subpsr = psr.add_subparsers(title="subcommands", dest='subcommand')
    psr_parent = ArgumentParser(add_help=False)

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


def create_iso_jobs(args):
    if args.sample=='cks':
        df = cksgaia.io.load_table('cks')
    else:
        print("Sample not defined: {}".format(args.sample))

    for i, row in df.iterrows():
        id_starname = row.id_starname
        outdir = "{}/{}".format(args.baseoutdir, id_starname)
        print "mkdir -p {}; run_cksgaia.py run-iso {} {} {} &> {}/run-iso.log".format(outdir, args.driver, id_starname, outdir, outdir)


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

