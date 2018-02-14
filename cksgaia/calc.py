import pandas as pd
import numpy as np
from collections import OrderedDict

import cksgaia.io
from cksgaia.errors import equad

SEED = 0  # reproducability
rand = np.random.RandomState(SEED)

albedo = 0.30


def update_planet_parameters(df):
    physmerge = df
    prad_iso = []
    prad_iso_err = []
    sinc_iso = []
    sinc_iso_err = []
    teq_iso = []
    teq_iso_err = []
    sma_iso = []
    sma_iso_err = []
    # flux_zams = []
    # flux_zams_err = []

    for i, row in physmerge.iterrows():

        # stellar radius
        unc_rad = np.mean([row['iso_srad_err1'], -row['iso_srad_err2']])
        # 8% floor on stellar radius uncertainty
        # if unc_rad / row['iso_srad'] < 0.08:
        #     unc_rad = 0.08 * row['iso_srad']
        e_rad = rand.normal(row['iso_srad'], unc_rad, size=10000)

        # stellar mass
        unc_mstar = np.mean([row['iso_smass_err1'], -row['iso_smass_err2']])
        # 5% floor on stellar mass uncertainty
        # if unc_mstar / row['iso_smass'] < 0.05:
        #     unc_mstar = 0.05 * row['iso_smass']
        e_mstar = rand.normal(row['iso_smass'], unc_mstar, size=10000)

        # period
        unc_per = np.mean([row['koi_period_err1'], -row['koi_period_err2']])
        if unc_per == 0.0:
            e_per = row['koi_period']
        else:
            e_per = rand.normal(row['koi_period'], unc_per, size=10000)

        # Teff
        unc_teff = np.mean([row['cks_steff_err1'], -row['cks_steff_err2']])
        e_teff = rand.normal(row['cks_steff'], unc_teff, size=10000)

        # Planet radius (earth radii)
        unc_ror = np.mean([row['koi_ror_err1'], -row['koi_ror_err2']])
        e_ror = rand.normal(row['koi_ror'], unc_ror, size=10000)

        unc_depth = np.mean([row['koi_depth'], row['koi_depth']])
        e_depth = rand.normal(row['koi_depth'], unc_depth, size=10000)

        prad = e_ror * e_rad * 109.245
        # prad = np.sqrt(depth*1e-6) * e_rad * 109.245
        prad_iso.append(np.median(prad))
        prad_iso_err.append(np.std(prad))

        a = (e_mstar * (e_per / 365.) ** 2.) ** (1. / 3.) / 0.00465
        teq = e_teff * np.sqrt(e_rad / (2.0 * a)) * (1.0 - albedo) ** (1. / 4.)
        teq_iso.append(np.median(teq))
        teq_iso_err.append(np.std(teq))

        # insolation flux
        a = (e_mstar * (e_per / 365.) ** 2.) ** (1. / 3.)
        sma_iso.append(np.median(a))
        sma_iso_err.append(np.std(a))

        sinc = (e_teff / 5778.) ** 4.0 * (e_rad / a) ** 2.0
        sinc_iso.append(np.median(sinc))
        sinc_iso_err.append(np.std(sinc))

    physmerge['iso_prad'] = prad_iso
    physmerge['iso_prad_err1'] = np.array(prad_iso_err)
    physmerge['iso_prad_err2'] = -np.array(prad_iso_err)

    physmerge['iso_insol'] = sinc_iso
    physmerge['iso_insol_err1'] = np.array(sinc_iso_err)
    physmerge['iso_insol_err2'] = -np.array(sinc_iso_err)

    physmerge['iso_teq'] = teq_iso
    physmerge['iso_teq_err1'] = np.array(teq_iso_err)
    physmerge['iso_teq_err2'] = -np.array(teq_iso_err)

    physmerge['iso_sma'] = sma_iso
    physmerge['iso_sma_err1'] = np.array(sma_iso_err)
    physmerge['iso_sma_err2'] = -np.array(sma_iso_err)

    return physmerge


def table_statistics():
    d = OrderedDict()
    cache = 1
    df = cksgaia.io.load_table('cks', cache=cache)
    d['nstars-cks'] = "{}".format(len(df))

    df = cksgaia.io.load_table('cks+nea', cache=cache)
    d['ncand-cks'] = "{}".format(len(df))

    # Model dep errors
    for k in 'smass srad-dw srad-gi sage'.split():
        d['equad-' + k] = "{:.0f}\%".format(equad[k] * 100)

    # Error summaries with and without floor
    for table in ['iso', 'iso-floor']:
        df = cksgaia.io.load_table(table, cache=cache)
        for k in 'smass srad slogage'.split():
            if k == 'slogage':
                err = df['iso_' + k + '_err1']
                fmt = ".2f"
                unit = "~dex"
            else:
                fmt = ".1f"
                err = df['iso_' + k + '_frac_err'] * 100
                unit = "\%"
            for p in [5, 50, 95]:
                fmt2 = "{:%s}%s" % (fmt, unit)
                percentile = fmt2.format(np.percentile(err, p))
                d["{}-{}-err-{:02d}".format(table, k, p, unit)] = percentile

    df = cksgaia.io.load_table('cks+nea+iso-floor', cache=cache)
    d['cks-rp-frac-err-median'] = "{:.0f}\%".format(
        (df['iso_prad_err1'] / df['iso_prad']).median() * 100)
    d['cks-sinc-frac-err-median'] = "{:.0f}\%".format(
        (df['iso_insol_err1'] / df['iso_insol']).median() * 100)
    d['cks-sma-frac-err-median'] = "{:.1f}\%".format(
        (df['iso_sma_err1'] / df['iso_sma']).median() * 100)

    for k, v in d.iteritems():
        print r"{{{}}}{{{}}}".format(k, v)
