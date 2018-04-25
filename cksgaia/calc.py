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
    prad_giso_err = []
    sinc_iso = []
    sinc_giso_err = []
    teq_iso = []
    teq_giso_err = []
    sma_iso = []
    sma_giso_err = []
    # flux_zams = []
    # flux_zams_err = []

    for i, row in physmerge.iterrows():

        # stellar radius
        unc_rad = np.mean([row['giso_srad_err1'], -row['giso_srad_err2']])
        # 8% floor on stellar radius uncertainty
        # if unc_rad / row['giso_srad'] < 0.08:
        #     unc_rad = 0.08 * row['giso_srad']
        e_rad = rand.normal(row['giso_srad'], unc_rad, size=10000)

        # stellar mass
        unc_mstar = np.mean([row['giso_smass_err1'], -row['giso_smass_err2']])
        # 5% floor on stellar mass uncertainty
        # if unc_mstar / row['giso_smass'] < 0.05:
        #     unc_mstar = 0.05 * row['giso_smass']
        e_mstar = rand.normal(row['giso_smass'], unc_mstar, size=10000)

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
        prad_giso_err.append(np.std(prad))

        a = (e_mstar * (e_per / 365.) ** 2.) ** (1. / 3.) / 0.00465
        teq = e_teff * np.sqrt(e_rad / (2.0 * a)) * (1.0 - albedo) ** (1. / 4.)
        teq_iso.append(np.median(teq))
        teq_giso_err.append(np.std(teq))

        # insolation flux
        a = (e_mstar * (e_per / 365.) ** 2.) ** (1. / 3.)
        sma_iso.append(np.median(a))
        sma_giso_err.append(np.std(a))

        sinc = (e_teff / 5778.) ** 4.0 * (e_rad / a) ** 2.0
        sinc_iso.append(np.median(sinc))
        sinc_giso_err.append(np.std(sinc))

    physmerge['giso_prad'] = prad_iso
    physmerge['giso_prad_err1'] = np.array(prad_giso_err)
    physmerge['giso_prad_err2'] = -np.array(prad_giso_err)

    physmerge['giso_insol'] = sinc_iso
    physmerge['giso_insol_err1'] = np.array(sinc_giso_err)
    physmerge['giso_insol_err2'] = -np.array(sinc_giso_err)

    physmerge['giso_teq'] = teq_iso
    physmerge['giso_teq_err1'] = np.array(teq_giso_err)
    physmerge['giso_teq_err2'] = -np.array(teq_giso_err)

    physmerge['giso_sma'] = sma_iso
    physmerge['giso_sma_err1'] = np.array(sma_giso_err)
    physmerge['giso_sma_err2'] = -np.array(sma_giso_err)

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
                err = df['giso_' + k + '_err1']
                fmt = ".2f"
                unit = "~dex"
            else:
                fmt = ".1f"
                err = df['giso_' + k + '_frac_err'] * 100
                unit = "\%"
            for p in [5, 50, 95]:
                fmt2 = "{:%s}%s" % (fmt, unit)
                percentile = fmt2.format(np.percentile(err, p))
                d["{}-{}-err-{:02d}".format(table, k, p, unit)] = percentile

    df = cksgaia.io.load_table('cks+nea+iso-floor', cache=cache)
    d['cks-rp-frac-err-median'] = "{:.0f}\%".format(
        (df['giso_prad_err1'] / df['giso_prad']).median() * 100)
    d['cks-sinc-frac-err-median'] = "{:.0f}\%".format(
        (df['giso_insol_err1'] / df['giso_insol']).median() * 100)
    d['cks-sma-frac-err-median'] = "{:.1f}\%".format(
        (df['giso_sma_err1'] / df['giso_sma']).median() * 100)

    for k, v in d.iteritems():
        print r"{{{}}}{{{}}}".format(k, v)
