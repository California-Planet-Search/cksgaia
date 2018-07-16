# CKS-Gaia

A re-analysis of the CKS dataset using Gaia constrained parallaxes.

## Installation

Download data tables

1. xmatch_m17_gaiadr2-result.csv # CKS-Gaia crossmatch
2. kic_q0_q17.hdf # Table of noise properties
3. kepler_stellar17.csv.gz # Mathur 2017 stellar properties tables

## Cross-match the CKS stars with Gaia IDs.

Follow instructions [here](docs/gaia-xmatch.md)

## Isochrone analysis

Creates `data/isoclassify-grid.csv` and `data/isoclassify-direct.csv` input files

```
bin/run_cksgaia.py create-iso-batch 
```

Then create batch files

```
isoclassify batch direct data/isoclassify-direct.csv -o isoclassify/direct/ > isoclassify-direct.tot
isoclassify batch grid data/isoclassify-grid.csv -o isoclassify/grid/ > isoclassify-grid.tot
```

Then run in parallel

```
# head isoclassify-direct.tot | parallel # useful for testing
# head isoclassify-grid.tot | parallel # useful for testing
parallel :::: isoclassify-direct.tot
parallel :::: isoclassify-grid.tot
```

Scrape through the output director to create stellar parameters table.

```
$ run_cksgaia.py create-iso-table isoclassify <outputdir> <outfile>
```
e.g.
```
$ run_cksgaia.py create-iso-table isoclassify isocla-j17-gaia isocla-j17-gaia.csv
```


# Updating plots tables and values in paper
```
$ run_cksgaia.py create-table all -d ./ # Make Tables
$ run_cksgaia.py create-val all -d ./   # Make values for latex
$ run_cksgaia.py create-plots all -d ./ # Make figures
```

###


```
run_cksgaia.py create-iso-batch 
isoclassify batch direct data/isoclassify-direct.csv -o isoclassify/direct/ > isoclassify-direct.tot
isoclassify batch grid data/isoclassify-grid.csv -o isoclassify/grid/ > isoclassify-grid.tot
mkdir -p isoclassify/direct//K00001;isoclassify run direct K00001 --outdir isoclassify/direct//K00001 --csv data/isoclassify-direct.csv 
head isoclassify-direct.tot | parallel # useful for testing
head isoclassify-grid.tot | parallel # useful for testing
parallel :::: isoclassify-direct.tot # Takes about XX min on Erik's laptop to run
parallel :::: isoclassify-grid.tot # Takes about XX min on Erik's laptop to run
```