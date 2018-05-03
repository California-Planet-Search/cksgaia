import pylab as pl
import matplotlib
import numpy as np
import pandas as pd

import cksgaia.io

np.random.seed(0)

gap_corners = [(1, 1.5), (100, 2.0)]
sn_corners = [(3, 2.0), (100, 4.0)]
e_corners = [(0.5, 0.7), (30, 1.5)]

wid_bad = 0.4
wid_best = 0.6
sn_cen = 2.4
e_cen = 1.2

def plot_box(corners, *args, **kwargs):
    pl.plot([corners[0][0], corners[1][0]],
            [corners[0][1], corners[0][1]], *args, **kwargs)
    pl.plot([corners[0][0], corners[0][0]],
            [corners[0][1], corners[1][1]], *args, **kwargs)
    pl.plot([corners[0][0], corners[1][0]],
            [corners[1][1], corners[1][1]], *args, **kwargs)
    pl.plot([corners[1][0], corners[1][0]],
            [corners[0][1], corners[1][1]], *args, **kwargs)


def plot_dist(real_sample):
    alphascale = 1/real_sample['giso_prad_err1']
    alphascale /= np.max(alphascale)
    real_sample['alpha'] = alphascale

    for i, row in real_sample.iterrows():
        p = row['koi_period']
        r = row['giso_prad']
        e = row['giso_prad_err1']
        a = row['alpha']
        pl.errorbar(p, r, yerr=e, alpha=a, fmt='k.')
    # pl.errorbar(real_sample['koi_period'], real_sample['giso_prad'],
    #             yerr=real_sample['giso_prad_err1'], fmt='k.')
    pl.loglog()
    cksgaia.plot.config.logscale_rad(axis='y')

    ax = pl.gca()
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%0.1f'))
    pl.xlim(0.3, 100)
    pl.xlabel('Orbital Period [days]')
    # pl.axhline(1.75, lw=3, color='r', linestyle='dashed')

    plot_box(e_corners, 'b-', lw=2)
    plot_box(sn_corners, 'b-', lw=2)
    plot_box(gap_corners, 'g-', lw=2)

    xt = [0.3, 1, 3, 10, 30, 100]
    yt = [0.7, 1, 1.5, 2, 3, 4, 6]
    pl.xticks(xt, xt)
    pl.yticks(yt, yt)

    pl.minorticks_off()
    pl.grid(False)


def planets_in_boxes(sample, box1, box2, box3):
    num_in_boxes = []
    for box in [box1, box2, box3]:
        radlim_low = box[0][1]
        radlim_high = box[1][1]
        plim_low = box[0][0]
        plim_high = box[1][0]
        q = sample.query('@radlim_low <= giso_prad < @radlim_high & @plim_low < koi_period <= @plim_high')
        n = len(q)

        num_in_boxes.append(n)

    return num_in_boxes


def make_mock_sample(real_sample, e_cen, sn_cen, wid):
    mock_sample = real_sample.copy()
    subnep = mock_sample.query('1.75 < giso_prad <= 4.0')
    earths = mock_sample.query('1.0 <= giso_prad <= 1.75')

    subnep['giso_prad'] = sn_cen + np.random.normal(loc=0, scale=subnep['giso_prad_err1'],
                                                    size=len(subnep)) \
                          + np.random.uniform(sn_cen * (-wid / 2), sn_cen * (wid / 2), size=len(subnep))

    earths['giso_prad'] = e_cen + np.random.normal(loc=0, scale=earths['giso_prad_err1'],
                                                   size=len(earths)) \
                          + np.random.uniform(e_cen * (-wid / 2), e_cen * (wid / 2), size=len(earths))

    mock_sample = pd.concat([subnep, earths])

    return mock_sample


def figure_of_merit(real_boxes, test_boxes):
    real_boxes = np.array(real_boxes)
    test_boxes = np.array(test_boxes)

    chi = np.sum((real_boxes - test_boxes) ** 2 / (np.sqrt(real_boxes) + np.sqrt(test_boxes)))

    return chi


def wid_sim_plot():
    real_sample = cksgaia.io.load_table('cksgaia-planets-weights')

    mock_bad = make_mock_sample(real_sample, e_cen, sn_cen, wid_bad)
    mock_good = make_mock_sample(real_sample, e_cen, sn_cen, wid_best)

    fig = pl.figure(figsize=(10, 3.5))
    pl.subplots_adjust(bottom=0.2, wspace=0.2, left=0.07,
                       right=0.95, top=0.9)

    pl.subplot(1, 3, 1)
    plot_dist(real_sample)
    pl.title("Real detections")
    pl.xlabel('')

    pl.subplot(1, 3, 2)
    plot_dist(mock_good)
    pl.title("Simulation (width={}%)".format(wid_best*100))
    ax = pl.gca()
    ax.yaxis.set_ticklabels([""])
    pl.ylabel('')

    pl.subplot(1, 3, 3)
    plot_dist(mock_bad)
    pl.title("Simulation (width={}%)".format(wid_bad*100))
    ax = pl.gca()
    ax.yaxis.set_ticklabels([""])
    pl.ylabel('')
    pl.xlabel('')
