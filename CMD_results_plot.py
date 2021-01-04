
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii


"""
After selecting a pyUPMASK method, find the three worst performers for each
set (PHOT, PM). This condition is estimated by subtracting the nine metrics
in pyUPMASK from UPMASK, and then adding the results. The three smallest values
point to the three worst performers.
"""

# Select a method here
method = 'Voron'
pyU_prob_min = .9
UP_prob_min = .9
print("Method: {}".format(method))

metrics_n = [
    'LSR', 'BSL', 'HMS', 'MCC_5', 'TPR_5', 'PPV_5', 'MCC_9', 'TPR_9', 'PPV_9']
outp_py, outp_UP = 'output/pyUPMASK_600/', 'output/UPMASK_600/'

for ftype in ('PHOT', 'PM'):
    mfile_UP = 'metrics/metrics_' + ftype + '_UPMASK_600_25_H_auto.dat'
    mfile_py = 'metrics/metrics_' + ftype + '_' + method + '_H_auto.dat'

    met_U = ascii.read(Path(mfile_UP))
    met_pyU = ascii.read(Path(mfile_py))
    met_deltas = []
    for met in metrics_n:
        met_deltas.append(met_pyU[met] - met_U[met])
    met_deltas = np.array(met_deltas)
    met_sum = met_deltas.sum(0)
    idx = np.argsort(met_sum)

    # Plot 3 worst files
    for j in range(3):

        print(met_deltas.T[idx[j]])
        if ftype == 'PHOT':
            worst_f = met_U['Name'][idx[j]].replace('_1_', '_1.').replace(
                '_2_', '_2.').replace('_3_', '_3.').replace('_4_', '_4.')
        else:
            aa = met_U['Name'][idx[j]].split('_')
            worst_f = aa[0] + '.' + aa[1] + '_' + aa[2] + '.' + aa[3] + '_' +\
                aa[4] + '.' + aa[5] + '_' + aa[6] + '.' + aa[7] + '_' +\
                aa[8] + '.' + aa[9]
        print(worst_f)

        fn_py = outp_py + method + '/' + ftype + '/' + worst_f + '.dat'
        fn_UP = outp_UP + ftype + '/' + worst_f + '.dat'
        data_pyupmask = ascii.read(fn_py)
        data_upmask = ascii.read(fn_UP)

        # Make plot
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

        plt.suptitle(worst_f + ' N={}'.format(msk_real.sum()))

        plt.subplot(231)
        msk_pyupmask = data_pyupmask['probs_final'] > pyU_prob_min
        plt.title('pyUPMASK, P>={}, N={}'.format(
            pyU_prob_min, msk_pyupmask.sum()))
        plt.scatter(
            data_pyupmask['x'], data_pyupmask['y'], zorder=-1, c='grey',
            marker='.', alpha=.3)
        plt.scatter(
            data_pyupmask['x'][msk_real], data_pyupmask['y'][msk_real],
            zorder=1, c='r', s=25)
        plt.scatter(
            data_pyupmask['x'][msk_pyupmask], data_pyupmask['y'][msk_pyupmask],
            s=15)

        plt.subplot(232)
        plt.scatter(
            data_pyupmask[BV_col], data_pyupmask['V'], zorder=-1, c='grey',
            marker='.', alpha=.3)
        plt.scatter(
            data_pyupmask[BV_col][msk_real], data_pyupmask['V'][msk_real],
            zorder=1, c='r', s=25)
        plt.scatter(
            data_pyupmask[BV_col][msk_pyupmask],
            data_pyupmask['V'][msk_pyupmask], s=15)
        plt.gca().invert_yaxis()

        plt.subplot(233)
        txt = (
            "LSR={:.2f}\nBSL={:.2f}\nHMS__={:.2f}"
            "\nMCC_5={:.2f}\nTPR_5={:.2f}\nPPV_5={:.2f}"
            "\nMCC_9={:.2f}\nTPR_9={:.2f}\nPPV_9={:.2f}").format(
                *met_pyU[
                    'LSR', 'BSL', 'HMS', 'MCC_5', 'TPR_5',
                    'PPV_5', 'MCC_9', 'TPR_9', 'PPV_9'][idx[j]])
        plt.text(-.1, .4, txt)
        plt.axis('off')

        plt.subplot(234)
        msk_upmask = data_upmask['probability'] >= UP_prob_min
        plt.title('UPMASK, P>={}, N={}'.format(UP_prob_min, msk_upmask.sum()))
        plt.scatter(
            data_upmask['x'], data_upmask['y'], zorder=-1, c='grey',
            marker='.', alpha=.3)
        plt.scatter(
            data_upmask['x'][msk_real], data_upmask['y'][msk_real],
            zorder=1, c='r', s=25)
        plt.scatter(
            data_upmask['x'][msk_upmask], data_upmask['y'][msk_upmask],
            s=15)

        plt.subplot(235)
        plt.scatter(
            data_upmask[BV_col], data_upmask['V'], zorder=-1, c='grey',
            marker='.', alpha=.3)
        plt.scatter(
            data_upmask[BV_col][msk_real], data_upmask['V'][msk_real],
            zorder=1, c='r', s=25)
        plt.scatter(
            data_upmask[BV_col][msk_upmask], data_upmask['V'][msk_upmask],
            s=15)
        plt.gca().invert_yaxis()

        plt.subplot(236)
        txt = (
            "LSR={:.2f}\nBSL={:.2f}\nHMS__={:.2f}"
            "\nMCC_5={:.2f}\nTPR_5={:.2f}\nPPV_5={:.2f}"
            "\nMCC_9={:.2f}\nTPR_9={:.2f}\nPPV_9={:.2f}").format(
                *met_U[
                    'LSR', 'BSL', 'HMS', 'MCC_5', 'TPR_5',
                    'PPV_5', 'MCC_9', 'TPR_9', 'PPV_9'][idx[j]])
        plt.text(-.1, .4, txt)
        plt.axis('off')

        plt.show()
