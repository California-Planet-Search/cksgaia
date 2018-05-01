from collections import OrderedDict
import cksgaia.io
import numpy as np

def val_stat(return_dict=False):
    d = OrderedDict()

    df = cksgaia.io.load_table('j17',cache=1)
    stars = df.groupby('id_kic',as_index=False).nth(0)
    cut = stars
    d['cks-star-count'] = len(stars)
    d['cks-mag-star-count'] = len(cut)

    df = cksgaia.io.load_table('m17+gaia2+j17',cache=1)
    stars = df.groupby('id_kic',as_index=False).nth(0)
    cut = stars
    d['cks-gaia-star-count'] = len(cut)
    d['cks-gaia-sparallax-ferr-median'] = "{:.1f}".format(cut.eval('gaia2_sparallax_err / gaia2_sparallax').median() * 100)

    d['cks-gaia-distmod-err-median'] = "{:.2f}".format(cut.eval('gaia2_sparallax_err / gaia2_sparallax').median())

    df = cksgaia.io.load_table('m17+gaia2+j17+ext',cache=0)
    d['ak-median'] = "{:.03f}".format(df.ext_ak.median())
    d['ak-min'] = "{:.03f}".format(df.ext_ak.min())
    d['ak-max'] = "{:.03f}".format(df.ext_ak.max())
    d['ak-median-err'] = "{:.03f}".format(df.ext_ak_err.median())
    d['ebv-median'] = "{:.03f}".format(df.ext_ebv.median())
    d['ebv-median-err'] = "{:.03f}".format(df.ext_ebv_err.median())
    d['mk-err-median'] = "{:.03f}".format(df.m17_kmag_err.median())
    
    # Properties of cks+gaia2 table
    df = cksgaia.io.load_table('cksgaia-planets',cache=1)
    cut = df
    ferr = cut.eval('0.5*(koi_ror_err1 - koi_ror_err2) / koi_ror')
    d['ror-ferr-median'] = "{:.1f}".format(100*ferr.median())
    ferr = cut.eval('0.5*(giso_srad_err1 - giso_srad_err2) / giso_srad')
    d['srad-ferr-median'] = "{:.1f}".format(100*ferr.median())
    ferr= cut.eval('0.5*(giso_prad_err1 - giso_prad_err2) / giso_prad')
    d['prad-ferr-median'] = "{:.1f}".format(100*ferr.median())


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
