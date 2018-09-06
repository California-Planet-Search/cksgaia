import pylab as pl
import numpy as np

import cksgaia.plot.sample

def apply_filters(physmerge, mkplot=False, verbose=False, textable=False):
    # slope and intercept for subgiant filter
    ls, li = 0.00025, 0.20

    def _bipanel(physmerge, nrow, ncol, plti, eloc=(12.0, 100), aloc=(0.95, 0.85), atxt='full CKS sample',
                 letters=True):

        if mkplot:
            left = plti
            right = plti + 1

            letters = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l']
            left_letter = letters[left - 1]
            right_letter = letters[right - 1]

            nea = physmerge.copy()
            nea['iso_prad'] = nea['koi_prad']
            nea['iso_prad_err1'] = nea['koi_prad_err1']

            pl.subplot(nrow, ncol, left)
            v_cks = cksgaia.plot.sample.simplehist(physmerge, annotate='%s (%d)' % (atxt, len(physmerge)),
                                               stacked=True, va_anno=False,
                                               eloc=eloc, aloc=aloc, fill_valley=False)
            if letters:
                pl.annotate(left_letter + ")", xy=(0.03, 0.9),
                            xycoords='axes fraction', weight='bold')
            pl.xlabel('')
            ymax_left = pl.ylim()[1]
            ax = pl.gca()
            plti += 1

            texline = "%s  &  %4.3f  \\\\" % (atxt, v_cks)
            if textable:
                f = open('tmp.tex', 'a')
                print >>f, texline
                f.close()

            print texline

            ylimits = [0, ymax_left]
            pl.subplot(nrow, ncol, left)
            pl.ylim(ylimits)

            return left + 1

    nrow = 11
    ncol = 1
    plti = 1

    if mkplot:
        figure = pl.figure(figsize=(8, 20))
        pl.subplots_adjust(hspace=0, wspace=0.25, top=0.98, bottom=0.05, right=0.96, left=0.15)

    if verbose:
        print "Initial catalogue = %d planets." % len(physmerge)

    plti = _bipanel(physmerge, nrow, ncol, plti)

    pre = len(physmerge)
    crop = physmerge[physmerge['cks_fp'] == False]
    post = len(crop)
    if verbose:
        print "False positive filter removes %d planets." % (pre - post)
    plti = _bipanel(crop, nrow, ncol, plti, atxt='false pos. removed')

    pre = len(crop)
    crop = crop[crop['kic_kepmag'] <= 14.2]
    post = len(crop)
    if verbose:
        print "Kp < 14.2 filter removes %d planets." % (pre - post)
    plti = _bipanel(crop, nrow, ncol, plti, eloc=(12.0, 70), atxt='$Kp < 14.2$')

    # pre = len(crop)
    # crop = crop[crop['gdir_srad'] / crop['gdir_srad_err1'] > 10]
    # post = len(crop)
    # if verbose:
    #     print "R$_{\star}$ > 10 $\sigma$" % (pre - post)
    # plti = _bipanel(crop, nrow, ncol, plti, eloc=(12.0, 70), atxt='$R_{\star} < 10 \sigma$')

    pre = len(crop)
    crop = crop[~(crop['fur17_rcorr_avg'] > 1.05)]
    post = len(crop)
    if verbose:
        print "Furlan+17 Rp correction < 5%% filter removes %d planets." % (pre - post)
    plti = _bipanel(crop, nrow, ncol, plti, eloc=(12.0, 70), atxt='dilution $\leq$ 5%% (%d)')

    # pre = len(crop)
    # crop = crop[crop['fur17_rcorr_avg'] <= 1.05]
    # if mkplot:
    #     pl.subplot(nrow,ncol,plti)
    #     v = cksrad.plotting.simplehist(crop, annotate='dilution $\leq$ 5%% (%d)' % len(crop), stacked=True, va_anno=True)
    # plti += 1
    # post = len(crop)
    # if verbose:
    #     print "Furlan+17 Rp correction < 5%% filter removes %d planets." % (pre-post)

    pre = len(crop)
    crop = crop[crop['gaia2_gflux_ratio'] < 1.1]
    post = len(crop)
    if verbose:
        print "GAIA dilution < 1.1 filter removes %d planets." % (pre - post)
    plti = _bipanel(crop, nrow, ncol, plti, eloc=(12.0, 60), atxt='$G_{\\rm blend} < 1.1$')


    pre = len(crop)
    crop = crop[crop['koi_impact'] <= 0.9]
    post = len(crop)
    if verbose:
        print "b < 0.9 filter removes %d planets." % (pre - post)
    plti = _bipanel(crop, nrow, ncol, plti, eloc=(12.0, 60), atxt='$b < 0.9$')

    # pre = len(crop)
    # crop = crop[crop['gdir_prad_err1']/crop['gdir_prad'] <= 0.1]
    # post = len(crop)
    # if verbose:
    #     print "R$_{p}$ precision < 10%" % (pre - post)
    # plti = _bipanel(crop, nrow, ncol, plti, eloc=(12.0, 60), atxt=r'$Rp_{\rm err} < 10%$')


    # pre = len(crop)
    # crop = crop[crop['vshaped'] >= 0.3]
    # pl.subplot(nrow,ncol,plti)
    # v = cksrad.plotting.simplehist(crop, annotate='v-shaped $\geq 0.3$ (%d)' % len(crop))
    # plti += 1
    # post = len(crop)
    # print "v-shaped >= 0.3 filter removes %d planets." % (pre-post)

    # pre = len(crop)
    # crop = crop[crop['gdir_prad_err1']/crop['gdir_prad'] <= 0.12]
    # if mkplot:
    #     pl.subplot(nrow,ncol,plti)
    #     v = cksrad.plotting.simplehist(crop, annotate='$\sigma(R_p)/R_p \leq 12$%% (%d)' % len(crop), stacked=True, va_anno=True)
    # plti += 1
    # post = len(crop)
    # if verbose:
    #     print "Rp_err < 12%% filter removes %d planets." % (pre-post)

    # pre = len(crop)
    # crop = crop[(crop['koi_period'] < 100)]  # & (crop['koi_period'] > 1)]
    # post = len(crop)
    # if verbose:
    #     print "Period filter removes %d planets." % (pre - post)
    # plti = _bipanel(crop, nrow, ncol, plti, eloc=(12.0, 60), atxt='$P < 100$ d')

    pre = len(crop)
    crop = crop[crop['gdir_srad'] <= 10 ** (ls * (crop['cks_steff'] - 5500) + li)]
    post = len(crop)
    if verbose:
        print "Subgiant filter removes %d planets." % (pre - post)
    plti = _bipanel(crop, nrow, ncol, plti, eloc=(12.0, 60), atxt='giant stars removed')

    pre = len(crop)
    crop = crop[(crop['cks_steff'] <= 6500) & (crop['cks_steff'] >= 4700)]
    post = len(crop)
    if verbose:
        print "Teff filter removes %d planets." % (pre - post)
    plti = _bipanel(crop, nrow, ncol, plti, eloc=(12.0, 60), aloc=(0.95, 0.85), atxt='4700 K < $T_{\\rm eff}$ < 6500 K')

    diff = crop['giso2_sparallax'] - crop['gaia2_sparallax']
    diff_err = np.sqrt(crop['gaia2_sparallax_err'] ** 2 + crop['giso2_sparallax_err1'] ** 2)
    pre = len(crop)
    crop = crop[(np.abs(diff / diff_err) < 4)]
    post = len(crop)
    if verbose:
        print "Parallax filter removes %d planets." % (pre - post)
    plti = _bipanel(crop, nrow, ncol, plti, eloc=(12.0, 60), aloc=(0.95, 0.85), atxt='$\pi_{\rm iso} \approx \pi_{\rm dir}')


    # pre = len(crop)
    # crop = crop[(crop['koi_snr'] >= 12)]
    # post = len(crop)
    # if verbose:
    #     print "SNR > 10 filter removes %d planets." % (pre-post)
    # plti = _bipanel(crop, nrow, ncol, plti, eloc=(8.0, 50), atxt='SNR > 10')

    if mkplot:
        # pl.subplot(nrow, ncol, plti-2)
        # pl.xticks([0.7, 1.0, 1.8, 3.5, 6.0, 10.0, 20.0])
        # pl.xlabel('Planet Size [Earth radii]', fontsize=28)

        pl.subplot(nrow, ncol, plti - 1)
        pl.xticks([0.7, 1.0, 1.8, 3.5, 6.0, 10.0, 20.0])
        pl.xlabel('Planet Size [Earth radii]', fontsize=28)

        # pl.annotate('Planet Size [Earth radii]', xy=(0.5, 0.03), xycoords='figure fraction',
        #                fontsize=28, horizontalalignment='center')
        pl.annotate('Number of Planets', xy=(0.04, 0.6), xycoords='figure fraction',
                    fontsize=28, horizontalalignment='center', rotation=90)
        # pl.annotate('Number of Planets', xy=(0.53, 0.6), xycoords='figure fraction',
        #                  fontsize=28, horizontalalignment='center', rotation=90)

    if verbose:
        print
        print "Final sample = %d planets." % len(crop)

    return crop
