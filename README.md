# cksgaia
Code for a reanalysis of CKS data with Gaia constraints

# Installation

Download a bunch of Milky Way dust models compiled by Bovy et al. (2015a)

https://github.com/jobovy/mwdust

and install

# Physical parameters

This part of the code follows the CKS-Physical codebase

## 1 Isochrone analysis

Create jobs that convert the spectroscopic parameters from Petigura et
al. (2017) to physical stellar parameters. Generate a list of shell scripts
for isoclassify.

```
python bin/run_cksgaia.py create-iso-jobs isoclassify cks isocla-j17-fakegaia | grep mkdir > isocla-j17-fakegaia.tot # for fake gaia data
```

Then run in parallel

```
$ parallel :::: isoclassify-cks.tot
```

Scrape through the output director to create stellar parameters table.

```
$ run_cksgaia.py create-iso-table <mode> <baseoutdir> <outfile>

e.g.

$ run_cksgaia.py create-iso-table isoclassify isocla-j17-fakegaia isocla-j17-fakegaia.csv

# This will put it in the right place for the table generation
$ run_cksgaia.py create-iso-table isochrones isocla+isochr-dsep data/isochrones.csv
```

## 2. Create Master CSV Table

Then create the overall table.

```
$ run_cksgaia.py create-merged-table
```

## 3. Update paper 

```
$ run_cksgaia.py tex-tables # LaTeX tables
$ run_cksgaia.py tex-stats # Summary statistics
$ run_cksgaia.py create-plots # Make figures
$ run_cksgaia.py update-paper # Move files to paper directory
```

###