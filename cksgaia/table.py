import os

import cksgaia.io
import cksgaia.fitting

def weight_table(lines='all'):
    physmerge = cksgaia.io.load_table('cksgaia-planets-weights')
    cols = [
        'id_koicand', 'koi_snr', 'det_prob',
        'tr_prob', 'weight'
    ]
    if lines == 'all':
        outstr = physmerge.to_latex(
            columns=cols, escape=False, header=False, index=False, 
            float_format='%4.2f'
        )
    else:
        outstr = physmerge.iloc[0:int(lines)].to_latex(
            columns=cols, escape=False, header=False, index=False, 
            float_format='%4.2f'
        )

    return outstr.split('\n')[2:-3]


def weight_table_machine():
    physmerge = cksgaia.io.load_table('cksgaia-planets-weights')

    full_cols = [
        'id_koicand', 'koi_snr', 'det_prob', 'tr_prob', 'weight'
    ]

    lines = []
    lines.append(", ".join(full_cols))

    for i, row in physmerge.iterrows():
        row_str = "{:s}, {:.2f}, {:.3f}, {:.4f}, {:.2f}".format(
            row['id_koicand'], row['koi_snr'], row['det_prob'],
            row['tr_prob'], row['weight']
        )
        lines.append(row_str)

    return lines


def bins_table():
    lines = []

    for i, rad in enumerate(cksgaia.fitting.Redges[:-1]):
        lines.append("%4.2f--%4.2f  &  %4.2f \\\\" % (rad, cksgaia.fitting.Redges[i + 1], cksgaia.fitting.efudge[i]))

    return lines


def filters_table():
    physmerge = cksgaia.io.load_table('fulton17')
    crop = cksgaia.io.apply_filters(physmerge, mkplot=True, textable=True)
    f = open('tmp.tex', 'r')
    lines = f.readlines()
    f.close()
    os.system('rm tmp.tex')
    lines = [l.replace('\n', '') for l in lines]
    return lines


def star():
    df = cksgaia.io.load_table('cksgaia-planets',cache=1)
    df = df.groupby('id_starname',as_index=False).nth(0)
    df = df.sort_values(by='id_starname')
    lines = []
    for i, row in df.iterrows():
        s = r""
        s+="{id_starname:s} & "
        s+="{cks_steff:0.0f} & "
        s+="{cks_smet:0.2f} & "
        s+="{m17_kmag:0.1f} & "
        s+="{gaia2_sparallax:0.2f} & "
        s+="{gdir_srad:0.2f} & "

        s+="{giso_smass:0.2f} & "
        s+="{giso_srad:0.2f} & "
        s+="{giso_srho:0.2f} & "
        s+="{giso_slogage:0.2f} & "

        s+="{giso2_sparallax:0.2f} & "
        s+=r"{gaia2_gflux_ratio:0.2f} & " 
        s+=r"{fur17_rcorr_avg:.3f} \\"
        s = s.format(**row)
        s = s.replace('nan','\\nodata')
        lines.append(s)

    return lines


def planet():
    df = cksgaia.io.load_table('cksgaia-planets',cache=1)
    df = df.sort_values(by='id_koicand')
    lines = []
    for i, row in df.iterrows():
        
        # Include errors
        '''
        s = r""
        s+=r"{id_koicand:s} & "
        s+=r"{koi_period:0.1f} & "
        s+=r"{koi_ror:.5f}_{{ {koi_ror_err2:.5f} }}^{{ +{koi_ror_err1:.5f} }} & "  
        s+=r"{gdir_prad:.2f}_{{ {gdir_prad_err2:.2f} }}^{{ +{gdir_prad_err1:.2f} }} & "  
        s+=r"{giso_sma:.5f}_{{ {giso_sma_err2:.5f} }}^{{ +{giso_sma_err1:.5f} }} & "  
        s+=r"{giso_insol:.0f}_{{ {giso_insol_err2:.0f} }}^{{ +{giso_insol_err1:.0f} }} \\  "
        '''


        s = r""
        s+=r"{id_koicand:s} & "
        s+=r"{koi_period:0.1f} & "
        s+=r"{koi_ror:.5f}  & "  
        s+=r"{gdir_prad:.2f} & "  
        s+=r"{giso_sma:.5f} & "  
        s+=r"{giso_insol:.0f} \\  "


        s = s.format(**row)
        s = s.replace('nan','\\nodata')
        lines.append(s)

    return lines


def star_machine():
    df = cksgaia.io.load_table('cksgaia-planets', cache=1)
    df = df.groupby('id_starname', as_index=False).nth(0)
    df = df.sort_values(by='id_starname')

    cols = ['id_starname',
            'cks_steff', 'cks_steff_err1', 'cks_steff_err2',
            'cks_smet', 'cks_smet_err1', 'cks_smet_err2',
            'm17_kmag', 'm17_kmag_err',
            'gaia2_sparallax', 'gaia2_sparallax_err',
            'gdir_srad', 'gdir_srad_err1', 'gdir_srad_err2',
            'giso_smass', 'giso_smass_err1', 'giso_smass_err2',
            'giso_slogage', 'giso_slogage_err1', 'giso_slogage_err2',
            'gaia2_gflux_ratio', 'fur17_rcorr_avg']

    df = df[cols]

    lines = []
    head = ",".join(cols)
    lines.append(head)
    for i, row in df.iterrows():
        l = "{id_starname:s},\
{cks_steff:.0f}, {cks_steff_err1:.0f}, {cks_steff_err2:.0f},\
{cks_smet:.2f},{cks_smet_err1:.2f},{cks_smet_err2:.2f},\
{m17_kmag:.3f},{m17_kmag_err:.3f},\
{gaia2_sparallax:.3f},{gaia2_sparallax_err:.3f},\
{gdir_srad:.3f},{gdir_srad_err1:.3f},{gdir_srad_err2:.3f},\
{giso_smass:.3f},{giso_smass_err1:.3f},{giso_smass_err2:.3f},\
{giso_slogage:.2f},{giso_slogage_err1:.2f},{giso_slogage_err2:.2f},\
{gaia2_gflux_ratio:.3f},{fur17_rcorr_avg:.4f}".format(**row)

        lines.append(l)

    return lines


def planet_machine():
    df = cksgaia.io.load_table('cksgaia-planets', cache=1)
    df = df.sort_values(by='id_koicand')

    cols = ['id_koicand',
            'koi_period', 'koi_period_err1', 'koi_period_err2',
            'koi_ror', 'koi_ror_err1', 'koi_ror_err2',
            'gdir_prad', 'gdir_prad_err1', 'gdir_prad_err2',
            'giso_sma', 'giso_sma_err1', 'giso_sma_err2',
            'giso_insol', 'giso_insol_err1', 'giso_insol_err2']

    lines = []
    head = ",".join(cols)
    lines.append(head)
    for i, row in df.iterrows():
        l = "{id_koicand:s},\
{koi_period:.9f}, {koi_period_err1:.9f}, {koi_period_err2:.9f},\
{koi_ror:.6f}, {koi_ror_err1:.6f}, {koi_ror_err2:.6f},\
{gdir_prad:.3f},{gdir_prad_err1:.3f},{gdir_prad_err2:.3f},\
{giso_sma:.5f},{giso_sma_err1:.5f},{giso_sma_err2:.5f},\
{giso_insol:.1f},{giso_insol_err1:.1f},{giso_insol_err2:.1f}".format(**row)

        lines.append(l)

    return lines
