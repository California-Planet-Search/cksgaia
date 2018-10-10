# CKS-Gaia

A re-analysis of the CKS dataset using Gaia constrained parallaxes.

## Installation

Download data tables

1. [xmatch_m17_gaiadr2-result.csv](https://www.dropbox.com/s/ctjmp1tfp3l6nuz/xmatch_m17_gaiadr2-result.csv?dl=0) # CKS-Gaia crossmatch
2. [kic_q0_q17.hdf](https://www.dropbox.com/s/54vz099idveo6mp/kic_q0_q17.hdf?dl=0) # Table of noise properties
3. [kepler_stellar17.csv.gz](https://www.dropbox.com/s/w8gwk0a2ieoju7b/kepler_stellar17.csv.gz?dl=0) # Mathur 2017 stellar properties tables

## Cross-match the CKS stars with Gaia IDs.

Follow instructions [here](docs/gaia-xmatch.md)

## Isochrone analysis

```
run_cksgaia.py create-iso-batch 
isoclassify batch direct data/isoclassify-direct.csv green18 -o isoclassify/direct/ > isoclassify-direct.tot
isoclassify batch grid data/isoclassify-grid-parallax-yes.csv green18 -o isoclassify/grid-parallax-yes/ > isoclassify-grid-parallax-yes.tot
isoclassify batch grid data/isoclassify-grid-parallax-no.csv green18 -o isoclassify/grid-parallax-no/ > isoclassify-grid-parallax-no.tot

# mkdir -p isoclassify/grid-parallax-yes//K00001;isoclassify run grid K00001 --outdir isoclassify/grid-parallax-yes//K00001 --csv data/isoclassify-grid-parallax-yes.csv --dust green18 &> isoclassify/grid-parallax-yes//K00001/output.log

head isoclassify-direct.tot | parallel # useful for testing
head isoclassify-grid-parallax-no.tot | parallel # useful for testing
head isoclassify-grid-parallax-yes.tot | parallel # useful for testing

cat isoclassify-direct.tot | parallel # useful for testing
cat isoclassify-grid-parallax-no.tot | parallel # takes about 20min on Erik's laptop
cat isoclassify-grid-parallax-yes.tot | parallel # takes about 20min on Erik's laptop

parallel :::: isoclassify-direct.tot # Takes about XX min on Erik's laptop to run
parallel :::: isoclassify-grid.tot # Takes about XX min on Erik's laptop to run

# check for missing stars
run_cksgaia.py create-iso-table

```

Sometimes there is an issue communicating with the Bayestar server executing the following code in serial quickly gathers the last problematic stars

```
for i in `grep "Response" isoclassify/direct/*/output.log | awk -F'/' '{print $3}' ` ;do eval `grep $i isoclassify-direct.tot` ;done 
for i in `grep "Response" isoclassify/grid-parallax-no/*/output.log | awk -F'/' '{print $3}' ` ;do eval `grep $i isoclassify-grid-parallax-no.tot` ;done
for i in `grep "Response" isoclassify/grid-parallax-yes/*/output.log | awk -F'/' '{print $3}' ` ;do eval `grep $i isoclassify-grid-parallax-yes.tot` ;done
```

## Plots, tables, and values in paper.
```
$ run_cksgaia.py create-table all -d ./ # Make Tables
$ run_cksgaia.py create-val all -d ./   # Make values for latex
$ run_cksgaia.py create-plots all -d ./ # Make figures
```

