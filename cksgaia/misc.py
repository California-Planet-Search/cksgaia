
import numpy as np


def frac_err(physmerge, x, col):
    """Fractional error

    Calculate median fractional uncertainty on a parameter in a DataFrame. Used to plot typical uncertainties'

    Args:
        physmerge (DataFrame): input catalogue
        x (float): x location where to calculate the relative uncertainty
        col (string): name of column

    Returns:
        tuple: (-err, +err)

    """

    frac_err = 1 + np.nanmedian(physmerge['%s_err1' % col]/physmerge[col])

    err1 = x - x / frac_err
    err2 = x * frac_err - x

    return (err1, err2)
