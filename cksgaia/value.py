from collections import OrderedDict
import cksgaia.io
import numpy as np

def val_stat(return_dict=False):
    d = OrderedDict()
    df = cksgaia.io.load_table('j17+m17+gaia2',cache=1)
    stars = df.groupby('id_kic',as_index=False).first()
    d['cks-star-count'] = len(stars)
    cut = stars.query('kic_kepmag < 14.2')
    d['cks-mag-star-count'] = len(cut)

    '''
    cut = stars.query('kic_kepmag < 14.2 & gaia2_angdist < 1 ')
    d['cks-gaia-star-count'] = len(cut)
    
    cut = stars.query('kic_kepmag < 14.2 & gaia2_angdist < 1 & gaia1_angdist_n1 < 1')
    d['cks-gaia-binary-count'] = len(cut.id_kic.drop_duplicates())
    '''

    '''
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

    '''
    # Comparison with Silva15
    comp = cksgaia.plot.compare.ComparisonRadius('srad-s15')
    d['cks-s15-count'] = comp.x3.count()
    d['cks-s15-srad-mean'] = comp.mean_string()
    d['cks-s15-srad-std'] = comp.std_string()

    # Comparison with Huber13
    comp = cksgaia.plot.compare.ComparisonRadius('srad-h13')
    d['cks-h13-count'] = comp.x3.count()
    d['cks-h13-srad-mean'] = comp.mean_string()
    d['cks-h13-srad-std'] = comp.std_string()

    # Comparison with Johnson17
    comp = cksgaia.plot.compare.ComparisonRadius('srad-j17')
    d['cks-j17-count'] = comp.x3.count()
    d['cks-j17-srad-mean'] = comp.mean_string()
    d['cks-j17-srad-std'] = comp.std_string()

    # Comparison with Gaia 
    comp = cksgaia.plot.compare.ComparisonRadius('srad-gaia2')
    d['cks-gaia2-count'] = comp.x3.count()
    d['cks-gaia2-srad-err-mean'] = "{:.1f}\%".format((comp.x2err[0]/  comp.x2).mean() * 100)
    d['cks-gaia2-srad-mean'] = comp.mean_string()
    d['cks-gaia2-srad-std'] = comp.std_string()


    lines = []
    for k, v in d.iteritems():
        line = r"{{{}}}{{{}}}".format(k,v)
        lines.append(line)

    if return_dict:
        return d

    return lines
