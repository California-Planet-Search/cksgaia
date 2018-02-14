import numpy as np
import astropysics.constants as C
import pandas as pd

import cksgaia.io
import cksgaia.fitting
from cksgaia.config import *


def fit_cdpp(kicselect):
    """Fit polynomial to CDPP values

    Args:
        kicselect (DataFrame): DataFrame of KIC catalogue. CDPP columns should be arrays
                               with value > 0 for quarters the star was observed.

    Returns:
        DataFrame : same as input with several columns added
    """

    kicselect['kic_nquarters'] = kicselect['CDPP3'].apply(lambda x: np.clip(np.sum(x > 0), 0, 17))
    kicselect['kic_tobs'] = kicselect['kic_nquarters'].values * 90

    kicselect['kic_cdpp3'] = kicselect['CDPP3'].apply(lambda x: np.median(x[x > 0]))
    kicselect['kic_cdpp6'] = kicselect['CDPP6'].apply(lambda x: np.median(x[x > 0]))
    kicselect['kic_cdpp12'] = kicselect['CDPP12'].apply(lambda x: np.median(x[x > 0]))

    # c3 = np.ma.masked_values(np.vstack(kicselect.as_matrix(columns=['CDPP3'])[:,0]), 0.0)
    # c6 = np.ma.masked_values(np.vstack(kicselect.as_matrix(columns=['CDPP6'])[:,0]), 0.0)
    # c12 = np.ma.masked_values(np.vstack(kicselect.as_matrix(columns=['CDPP12'])[:,0]), 0.0)
    # kicselect['kic_cdpp3'] = np.median(c3, axis=1)
    # kicselect['kic_cdpp6'] = np.median(c6, axis=1)
    # kicselect['kic_cdpp12'] = np.median(c12, axis=1)

    cdpp_arr = kicselect.as_matrix(columns=['kic_cdpp3', 'kic_cdpp6', 'kic_cdpp12']).transpose()
    pfit = np.polyfit(1 / np.sqrt([3., 6., 12.]), cdpp_arr, 2)

    kicselect['kic_cdpp_fit0'] = pfit[2, :]
    kicselect['kic_cdpp_fit1'] = pfit[1, :]
    kicselect['kic_cdpp_fit2'] = pfit[0, :]

    return kicselect


def detection_prob(prad, per, kicselect, nkic=None, step=False):
    if nkic is None:
        nkic = kicselect['id_kic'].count()

    # Put this planet around all other stars
    # Calculate new durations
    a = (C.G * kicselect['kic_smass'].values * C.Ms * ((per * (24 * 3600.)) / (2 * np.pi)) ** 2) ** (1 / 3.) * C.aupercm
    R = kicselect['kic_srad'].values * C.Rs * C.aupercm
    durations = (per * 24. / np.pi) * np.arcsin(R / a)
    rors = (prad * (C.Re / C.Rs)) / kicselect['kic_srad'].values

    x = 1 / np.sqrt(durations)
    cdpp_durs = kicselect['kic_cdpp_fit0'].values + kicselect['kic_cdpp_fit1'].values * x + kicselect[
                                                                                                'kic_cdpp_fit2'].values * x ** 2

    # Calculate SNR for other stars
    other_snr = rors ** 2 * (per / kicselect['kic_tobs'].values) ** -0.5 * (1 / (cdpp_durs * 1e-6))

    # s = cksrad.fitting.logistic(other_snr, step=step)
    s = cksgaia.fitting.gamma_complete(other_snr, step=step)
    s[np.isnan(s)] = 0.0
    det = np.sum(s) / np.float(nkic)

    return det


def get_weights(kois, kicselect):
    """
    Add completeness columns det_prob, tr_prob, and weight into
    kois DataFrame

    Args:
        kois (DataFrame): input dataframe with candidates
        kicselect (DataFrame): DataFrame with all KIC stars in sample

    Returns:
        DataFrame: same as kois with det_prob, tr_prob, and weight columns
                   added
    """

    x = 1 / np.sqrt(kois['koi_duration'].values)
    kois['koi_cdpp_dur'] = kois['kic_cdpp_fit0'].values + \
                           kois['kic_cdpp_fit1'].values * x + \
                           kois['kic_cdpp_fit2'].values * x ** 2
    kois['koi_snr'] = kois['koi_ror'].values ** 2 * \
                      (kois['koi_period'].values / kois['kic_tobs'].values) ** -0.5 * \
                      1 / (kois['koi_cdpp_dur'].values * 1e-6)

    det_prob = []
    jobs = []
    nkic = kicselect['id_kic'].count()

    pers = kois['koi_period'].values
    prads = kois['iso_prad'].values

    for i in range(len(pers)):
        per = pers[i]
        prad = prads[i]

        det = detection_prob(prad, per, kicselect, nkic=nkic)
        det_prob.append(det)

        # print i, per, prad, det

    det_prob = np.array(det_prob)
    # tr_prob = 0.7/kois['koi_dor'].values

    mstar_g = kois['iso_smass'].values * C.Ms
    per_s = kois['koi_period'].values * (24 * 3600.)

    a_cm = (C.G * mstar_g * (per_s / (2 * np.pi)) ** 2) ** (1 / 3.)
    # a = (num/(2*np.pi)**2)**(1/3.) * C.aupercm
    R_cm = kois['iso_srad'].values * C.Rs
    tr_prob = 0.7 * R_cm / a_cm

    weights = 1 / (det_prob * tr_prob)

    kois['det_prob'] = det_prob
    kois['tr_prob'] = tr_prob
    kois['weight'] = weights

    return kois


def weight_merge(physmerge):
    kic = cksgaia.io.load_table('kic')
    kicselect = cksgaia.io.load_table('kic-filtered')
    kicselect = fit_cdpp(kicselect)

    physmerge = pd.merge(physmerge, kic, on='id_kic')

    physmerge = fit_cdpp(physmerge)
    physmerge = get_weights(physmerge, kicselect)

    return physmerge