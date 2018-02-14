import numpy as np
import cksgaia.io
import pandas as pd

# Add radii in quadrature
equad = {
    'smass': 0.02,  # Based on comparison of results returned by Dartmouth models
    'srad-dw': 0.02,  # Same as mass. Supported by Huber+17 AS-TGAS dist. comparison
    'srad-gi': 0.10,  # Same as mass. Supported by Huber+17 AS-TGAS dist. comparison
    'sage': 0.10  # model-dependent error floor
}


def frac_equad(val, err1, err2, equad):
    """Add fractional errors in quadrature

    Args:
        val (float): value
        err1 (float): upper uncert
        err2 (float): lower uncert (negative number)
        efloor (float): Add fractional errors quadrature to be added in quad

    Returns:
        val, err1, err2

    """
    err1 = np.sqrt(err1 ** 2 + (val * equad) ** 2)  # Add in quad
    err2 = -1.0 * np.sqrt(err2 ** 2 + (val * equad) ** 2)  # Add in quad
    return val, err1, err2


def frac_efloor(val, err1, err2, efloor):
    """Compute errors after adopting an error floor

    Args:
        val (float): value
        err1 (float): upper uncert
        err2 (float): lower uncert (negative number)
        efloor (float): do not return errors lower than this
            fractional precision

    Returns:
        val, err1, err2

    """
    frac_err = _frac_err(val, err1, err2)
    frac_err = np.vstack([frac_err, np.ones_like(frac_err) * efloor])
    frac_err = np.max(frac_err, axis=0)
    err1 = frac_err * val
    err2 = -1.0 * val * frac_err / (1 + frac_err)
    return val, err1, err2


def frac_err(val, err1, err2):
    """Compute fractional errors given double-sided error bars"""
    frac_err1 = err1 / val
    frac_err2 = err2 / (val + err2)
    frac_err = np.array(0.5 * (frac_err1 - frac_err2))
    return frac_err


def add_equad(df):
    df.index = df.id_starname
    for key in ['smass', 'srad-dw', 'srad-gi', 'sage']:
        _equad = equad[key]
        if key == 'srad-dw':
            idx = df.query('iso_slogg > 3.9').index
            key = 'srad'
        elif key == 'srad-gi':
            idx = df.query('iso_slogg < 3.9').index
            key = 'srad'
        else:
            idx = df.index

        df2 = df.copy()
        df2 = df2.ix[idx]

        val = df2['iso_' + key]
        err1 = df2['iso_' + key + '_err1']
        err2 = df2['iso_' + key + '_err2']
        _, err1, err2 = cksgaia.errors.frac_equad(val, err1, err2, _equad)
        df.loc[idx, 'iso_' + key + '_err1'] = np.array(err1)
        df.loc[idx, 'iso_' + key + '_err2'] = np.array(err2)

    # Must update slogage
    df['iso_slogage_err1'] = (
        np.log10(df.iso_sage + df.iso_sage_err1) + 9 - df.iso_slogage
    )
    df['iso_slogage_err2'] = (
        np.log10(df.iso_sage + df.iso_sage_err2) + 9 - df.iso_slogage
    )
    return df


def add_frac_err(df):
    for key in ['srad', 'smass']:
        val = df['iso_' + key]
        err1 = df['iso_' + key + '_err1']
        err2 = df['iso_' + key + '_err2']
        df['iso_' + key + '_frac_err'] = frac_err(val, err1, err2)

    return df


def rchisq(df, key, prefixes):
    err1 = df['iso_' + key + '_err1']
    err2 = df['iso_' + key + '_err2']
    df['iso_' + key + '_err'] = 0.5 * (err1 - err2)
    pre1 = prefixes[0]
    pre2 = prefixes[1]
    x1 = df[pre1 + key]
    x1err = df[pre1 + key + '_err']
    x2 = df[pre2 + key]
    x2err = df[pre2 + key + '_err']
    df['norm_sq_diff'] = (x2 - x1) ** 2 / (x2err ** 2 + x1err ** 2)
    _rchisq = df['norm_sq_diff'].sum() / len(df)
    return _rchisq, df


