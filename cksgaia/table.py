import os

import cksgaia.io
import cksgaia.fitting

def weight_table(lines='all'):
    physmerge = cksgaia.io.load_table('cksgaia-planets-weights')
    cols = [
        'id_koicand', 'koi_period', 'iso_prad', 'koi_snr', 'det_prob', 
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

    return outstr.split('\n')

def weight_table_machine():
    physmerge = cksgaia.io.load_table('fulton17-weights')

    full_cols = [
        'id_koicand', 'koi_period', 'koi_period_err1', 'koi_period_err2',
        'iso_prad', 'iso_prad_err1', 'iso_prad_err2'
        'koi_snr', 'det_prob', 'tr_prob', 'weight'
    ]

    lines = []
    lines.append(", ".join(full_cols))

    for i, row in physmerge.iterrows():
        row_str = "{:s}, {:10.8f}, {:.1e}, {:.1e}, {:.2f}, {:.2f}, {:.2f}, {:.2f}, {:.3f}, {:.4f}, {:.2f}".format(
            row['id_koicand'], row['koi_period'], row['koi_period_err1'], row['koi_period_err2'],
            row['iso_prad'], row['iso_prad_err1'], row['iso_prad_err2'],
            row['koi_snr'], row['det_prob'], row['tr_prob'], row['weight']
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
        s+="{ext_ak:0.3f} & "
        s+="{giso_srad:0.2f} & "
        s+="{giso_smass:0.2f} & "
        s+="{giso_slogage:0.2f} & "
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
        s+=r"{giso_prad:.2f}_{{ {giso_prad_err2:.2f} }}^{{ +{giso_prad_err1:.2f} }} & "  
        s+=r"{giso_sma:.5f}_{{ {giso_sma_err2:.5f} }}^{{ +{giso_sma_err1:.5f} }} & "  
        s+=r"{giso_insol:.0f}_{{ {giso_insol_err2:.0f} }}^{{ +{giso_insol_err1:.0f} }} \\  "
        '''


        s = r""
        s+=r"{id_koicand:s} & "
        s+=r"{koi_period:0.1f} & "
        s+=r"{koi_ror:.5f}  & "  
        s+=r"{giso_prad:.2f} & "  
        s+=r"{giso_sma:.5f} & "  
        s+=r"{giso_insol:.0f} \\  "


        s = s.format(**row)
        s = s.replace('nan','\\nodata')
        lines.append(s)
    return lines
