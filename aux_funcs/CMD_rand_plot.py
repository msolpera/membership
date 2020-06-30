
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import ascii


"""
Picks a random file (output from the pyUPMASk and UPMASK analysis) and
generates a coordinates and CMD plot.

A method for pyUPMASK and for the H measure metric need to be selected.
"""

method = 'autoperc_inner_GUMM3'
Hmethod = 'auto'

print("Method: {}".format(method))

while True:
    # Select random file from pyUPMASK results folder
    feature = ('PHOT', 'PM')
    feature_ran = feature[np.random.choice(2)]
    fold = '../output/' + method + '/' + feature_ran + '/'

    files = os.listdir(fold)
    file = np.random.choice(files)
    print(file)

    # Read pyUPMASK data
    data_pyupmask = ascii.read(fold + file)

    # Load the same file from the UPMASK run
    UP_fold = '../TEST_SYNTH_CLUSTS/test_results/100_UPMASK_res/'
    data_upmask = ascii.read(UP_fold + file)

    # Load pyUPMASK metrics
    metrics_file = '../TEST_SYNTH_CLUSTS/test_results/metrics_' +\
        feature_ran + '_' + method + '_H_' + Hmethod + '.dat'
    metrics_data = ascii.read(metrics_file)
    msk = metrics_data['Name'] == file[:-4].replace('.', '_')
    met_pyU = metrics_data[msk]

    # Load UPMASK metrics
    metrics_file = '../TEST_SYNTH_CLUSTS/test_results/metrics_UP-' +\
        feature_ran + '_H_' + Hmethod + '.dat'
    metrics_data = ascii.read(metrics_file)
    msk = metrics_data['Name'] == file[:-4].replace('.', '_')
    met_U = metrics_data[msk]

    try:
        msk_real = data_upmask['field'] == 0
        BV_col = 'B_V'
    except KeyError:
        BV_col = 'BV'
        msk_real = []
        for _id in data_upmask['ID']:
            if str(_id)[0] == '1':
                msk_real.append(True)
            else:
                msk_real.append(False)
        msk_real = np.array(msk_real)

    plt.suptitle(file + ' N={}'.format(msk_real.sum()))

    plt.subplot(231)
    msk_pyupmask = data_pyupmask['probs_final'] > .0
    plt.title('pyUPMASK N={}'.format(msk_pyupmask.sum()))
    plt.scatter(
        data_upmask['x'][msk_real], data_upmask['y'][msk_real],
        zorder=1, c='r', s=25)
    plt.scatter(
        data_pyupmask['x'][msk_pyupmask], data_pyupmask['y'][msk_pyupmask],
        s=15)

    plt.subplot(232)
    plt.scatter(
        data_upmask[BV_col][msk_real], data_upmask['V'][msk_real],
        zorder=1, c='r', s=25)
    plt.scatter(
        data_pyupmask[BV_col][msk_pyupmask], data_pyupmask['V'][msk_pyupmask],
        s=15)
    plt.gca().invert_yaxis()

    plt.subplot(233)
    txt = (
        "LSR___={:.2f}  BSL___={:.2f}  HMS__={:.2f}"
        "\nMCC_5={:.2f}  TPR_5={:.2f}  PPV_5={:.2f}"
        "\nMCC_9={:.2f}  TPR_9={:.2f}  PPV_9={:.2f}").format(
            *np.array(met_pyU[
                'LSR', 'BSL', 'HMS', 'MCC_5', 'TPR_5',
                'PPV_5', 'MCC_9', 'TPR_9', 'PPV_9'])[0])
    plt.text(-.1, .4, txt, fontsize=15)
    plt.axis('off')

    plt.subplot(234)
    msk_upmask = data_upmask['probability'] >= .9
    plt.title('UPMASK, P>=.9, N={}'.format(msk_upmask.sum()))
    plt.scatter(
        data_upmask['x'][msk_real], data_upmask['y'][msk_real],
        zorder=1, c='r', s=25)
    plt.scatter(
        data_upmask['x'][msk_upmask], data_upmask['y'][msk_upmask],
        s=15)

    plt.subplot(235)
    plt.scatter(
        data_upmask[BV_col][msk_real], data_upmask['V'][msk_real],
        zorder=1, c='r', s=25)
    plt.scatter(
        data_upmask[BV_col][msk_upmask], data_upmask['V'][msk_upmask], s=15)
    plt.gca().invert_yaxis()

    plt.subplot(236)
    txt = (
        "LSR___={:.2f}  BSL___={:.2f}  HMS__={:.2f}"
        "\nMCC_5={:.2f}  TPR_5={:.2f}  PPV_5={:.2f}"
        "\nMCC_9={:.2f}  TPR_9={:.2f}  PPV_9={:.2f}").format(
            *np.array(met_U[
                'LSR', 'BSL', 'HMS', 'MCC_5', 'TPR_5',
                'PPV_5', 'MCC_9', 'TPR_9', 'PPV_9'])[0])
    plt.text(-.1, .4, txt, fontsize=15)
    plt.axis('off')

    plt.show()
