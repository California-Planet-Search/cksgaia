from collections import OrderedDict
import cksgaia.io

# Code from EAP's K2-24 paper use as a template
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
    lines = []
    for k, v in d.iteritems():
        line = r"{{{}}}{{{}}}".format(k,v)
        lines.append(line)

    if return_dict:
        return d

    return lines
