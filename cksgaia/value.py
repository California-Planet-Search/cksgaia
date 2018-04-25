from collections import OrderedDict
import cksgaia.io
import numpy as np

def val_stat(return_dict=False):
    d = OrderedDict()
    df = cksgaia.io.load_table('j17+m17+gaia1',cache=1)
    stars = df.groupby('id_kic',as_index=False).first()
    d['cks-star-count'] = len(stars)
    cut = stars.query('kic_kepmag < 14.2')
    d['cks-mag-star-count'] = len(cut)

    cut = stars.query('kic_kepmag < 14.2 & gaia1_angdist < 1 ')
    d['cks-gaia-star-count'] = len(cut)
    
    cut = stars.query('kic_kepmag < 14.2 & gaia1_angdist < 1 & gaia1_angdist_n1 < 1')
    d['cks-gaia-binary-count'] = len(cut.id_kic.drop_duplicates())

    df = cksgaia.io.load_table('j17+m17+extinct',cache=0)
    df = df.query('kic_kepmag < 14.2')
    d['ak-median'] = "{:.03f}".format(df.ak.median())
    d['ak-median-err'] = "{:.03f}".format(df.ak_err.median())
    d['ebv-median'] = "{:.03f}".format(df.ebv.median())
    d['ebv-median-err'] = "{:.03f}".format(df.ebv_err.median())
    d['mk-err-median'] = "{:.03f}".format(df.m17_kmag_err.median())

    # Properties of cks+gaia2 table
    df = cksgaia.io.load_table('fakegaia-merged',cache=1)
    cut = df.query('kic_kepmag < 14.2')
    ferr= cut.eval('0.5*(koi_ror_err1 - koi_ror_err2) / koi_ror')
    d['ror-ferr-median'] = "{:.1f}".format(100*np.nanmedian(ferr))
    ferr = cut.eval('0.5*(iso_srad_err1 - iso_srad_err2) / iso_srad')
    d['srad-ferr-median'] = "{:.1f}".format(100*np.nanmedian(ferr))
    ferr= cut.eval('0.5*(iso_prad_err1 - iso_prad_err2) / iso_prad')
    d['prad-ferr-median'] = "{:.1f}".format(100*np.nanmedian(ferr))

    # Comparison with AS
    df = cksgaia.io.load_table('cks+gaia2+s15',cache=1)
    logdiff = np.log10(df.iso_srad) - np.log10(df.s15_srad)
    d['cks-s15-count'] = "{}".format(logdiff.count())
    d['cks-s15-srad-mean'] = "{:.1f}".format(100*(10**logdiff.mean() - 1))
    d['cks-s15-srad-std'] = "{:.1f}".format(100*(10**logdiff.std() - 1))
    ferr = df.eval('0.5 * (s15_srad_err1 - s15_srad_err2) / s15_srad')
    d['s15-srad-ferr-median'] = "{:.1f}".format(100*ferr.median())

    df = cksgaia.io.load_table('cks+gaia2+h13',cache=1)
    logdiff = np.log10(df.iso_srad) - np.log10(df.h13_srad) 
    d['cks-h13-count'] = "{}".format(logdiff.count())
    d['cks-h13-srad-mean'] = "{:.1f}".format(100*(10**logdiff.mean() - 1))
    d['cks-h13-srad-std'] = "{:.1f}".format(100*(10**logdiff.std() - 1))
    ferr = df.eval('h13_srad_err / h13_srad')
    d['h13-srad-ferr-median'] = "{:.1f}".format(100*ferr.median())

    lines = []
    for k, v in d.iteritems():
        line = r"{{{}}}{{{}}}".format(k,v)
        lines.append(line)

    if return_dict:
        return d

    return lines
