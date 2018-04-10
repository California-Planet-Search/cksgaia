# Cross-Matching Gaia

We have two methods for cross-matching Gaia catalogs


## Method 1: ADQL at Gaia Archive

Validated against Gaia DR1. For DR2 change gaiadr1 -> gaiadr2

## CKS

1. Generate VOTables and upload to Gaia Archive
2. Log into Gaia archive as `epetigur`
3. Perform the following joins that query all gaia sources within 8 arcsec of each CKS source, or stellar17 source.
4. Give the job a sensible name like `xmatch_cks_gaiadr1`
5. Run the query
6. Download results as a votable
7. Move into data

```
SELECT *,distance(
  POINT('ICRS', cks.m17_ra, cks.m17_dec),
  POINT('ICRS', gaia.ra, gaia.dec)) AS dist
FROM gaiadr1.gaia_source AS gaia, user_epetigur.cks AS cks
WHERE 1=CONTAINS(
  POINT('ICRS', cks.m17_ra, cks.m17_dec),
  CIRCLE('ICRS', gaia.ra, gaia.dec, 0.00222)
)
```

## Mathur17

Same as CKS, but with this query

```
SELECT *,distance(
  POINT('ICRS', m17.m17_ra, m17.m17_dec),
  POINT('ICRS', gaia.ra, gaia.dec)) AS dist
FROM gaiadr1.gaia_source AS gaia, user_epetigur.m17 AS m17
WHERE 1=CONTAINS(
  POINT('ICRS', m17.m17_ra, m17.m17_dec),
  CIRCLE('ICRS', gaia.ra, gaia.dec, 0.00222)
)
```

## Method 2: XMatch at CDS

upload to xmatch. cross-reference with Gaia (search size of 8 arcsec)

save as

data/cks-xmatch-results.csv




