# CKS-Gaia

A re-analysis of the CKS dataset using Gaia constrained parallaxes.

## Installation

Download data tables

1. xmatch_m17_gaiadr2-result.csv # CKS-Gaia crossmatch
2. kic_q0_q17.hdf # Table of noise properties
3. kepler_stellar17.csv.gz # Mathur 2017 stellar properties tables


## Cross-match the CKS stars with Gaia IDs.

Follow instructions [here](docs/gaia-xmatch.md)

## Physical parameters

This part of the code follows the CKS-Physical codebase

## 1 Isochrone analysis

Create jobs that convert the spectroscopic parameters from Petigura et
al. (2017) to physical stellar parameters. Generate a list of shell scripts
for isoclassify.

```
python bin/run_cksgaia.py create-iso-jobs isoclassify cks+gaia2 <outputdir> | grep mkdir > <outputdir>/isoclassify-cksgaia.tot
```

Then run in parallel

```
$ cd <outputdir>
$ parallel :::: isoclassify-cks.tot
```

Scrape through the output director to create stellar parameters table.

```
$ run_cksgaia.py create-iso-table isoclassify <outputdir> <outfile>

e.g.

$ run_cksgaia.py create-iso-table isoclassify isocla-j17-gaia isocla-j17-gaia.csv



# Updating plots tables and values in paper
```
$ run_cksgaia.py create-table all -d ./ # Make Tables
$ run_cksgaia.py create-val all -d ./   # Make values for latex
$ run_cksgaia.py create-plots all -d ./ # Make figures
```
