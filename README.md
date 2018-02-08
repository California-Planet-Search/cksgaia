# cksgaia
Code for a reanalysis of CKS data with Gaia constraints

# Installation

Download a bunch of Milky Way dust models compiled by Bovy et al. (2015a)

https://github.com/jobovy/mwdust

and install


# Physical parameters

This part of the code follows the CKS-Physical codebase

## 1 Isochrone analysis

Create jobs that convert the spectroscopic parmaeters stored in
cks_v2.csv into physical parameters. Generate a list of shell scripts
for running morton's isochrones code. On the isoclassify branch we can
run Huber's isoclassify code

python bin/run_cksgaia.py create-iso-jobs isoclassify cks | grep mkdir > isoclassify-cks.tot

```
$ parallel :::: isoclassify-cks.tot
```

Scrape through the isochrone output files to create the isochrone table. 

```
$ run_cksgaia.py create-iso-table <mode> <baseoutdir> <outfile>
e.g.

$ run_cksgaia.py create-iso-table isochrones isocla+isochr-mist isocla+isochr-mist.csv

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