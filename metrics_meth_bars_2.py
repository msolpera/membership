
# from pathlib import Path
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from auxFuncs import tie_min, tie_max, WinTieLoss


# Folder where the files are located
fold = "metrics/"


def main():
    """
    Plot the vertical bar plots for all the statistics, separated into PM,
    PHOT, and combined results.
    """
    configs = ("Voron", "kNNde", "Agglo", "MiniB", "KMean", "Gauss")
    Hval = 'auto'  # 'symm', 'SR05'
    N_UPMASK = "25"  # "50"
    ASteCA = 'true'  # 'true'

    results = []
    # results_phot = []
    for i, m in enumerate(configs):
        print(m)
        if ASteCA == 'true':
            pyUP_PHOT, pyUP_PM, AS_PHOT, AS_PM = readTables(
                N_UPMASK, Hval, m, flagasteca=True)
            win_PHOT, loss_PHOT, emp_PHOT, win_PM, loss_PM, emp_PM =\
                WinTieLoss(pyUP_PHOT, pyUP_PM, AS_PHOT,
                           AS_PM, 'metrics_bars')
            results.append(makePlot(
            N_UPMASK, tie_max, tie_min, m, win_PHOT, loss_PHOT,
            emp_PHOT, win_PM, loss_PM, emp_PM, flagasteca=True))
        else:
            pyUP_PHOT, pyUP_PM, UP_PHOT, UP_PM = readTables(
                N_UPMASK, Hval, m)
            win_PHOT, loss_PHOT, emp_PHOT, win_PM, loss_PM, emp_PM =\
                WinTieLoss(pyUP_PHOT, pyUP_PM, UP_PHOT,
                           UP_PM, 'metrics_bars')
            results.append(makePlot(
                N_UPMASK, tie_max, tie_min, m, win_PHOT, loss_PHOT,
                emp_PHOT, win_PM, loss_PM, emp_PM))

    fig = plt.figure(figsize=(22, 8))
    conf = ("VOR", "KNN", "AGG", "MBK", "KMS", "GMM")
    for j, l in enumerate(conf):
        ax = plt.subplot(2, 3, j+1)
        barsPlot(ax, results[j][0], results[j][1])
        if j == 0:
            plt.title(l, loc='right', fontsize=12)
        else:
            plt.title(l, fontsize=12)
        ax.invert_xaxis()
        if j != 0:
            ax.legend().remove()
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
    fig.tight_layout()
    plt.savefig('metrics_bars_PM.png', dpi=300, bbox_inches='tight')

    # plt.show()


def readTables(N_UPMASK, H, m, flagasteca=False):
    """
    """
    pyUP_PHOT = Table.read(
        fold + 'metrics_PHOT_{}_H_{}.dat'.format(m, H), format='ascii')
    pyUP_PM = Table.read(
        fold + 'metrics_PM_{}_H_{}.dat'.format(m, H), format='ascii')

    UP_PHOT = Table.read(
        fold + '/metrics_PHOT_UP_600_' + str(N_UPMASK) + '_H_' + H +
        '.dat', format='ascii')
    UP_PM = Table.read(
        fold + '/metrics_PM_UP_600_' + str(N_UPMASK) + '_H_' + H +
        '.dat', format='ascii')
    if flagasteca == True:
        AS_PHOT = Table.read(
            fold + '/metrics_PHOT_ASteCA' + '_H_' + H +
            '.dat', format='ascii')

        AS_PM = Table.read(
            fold + '/metrics_PM_ASteCA' + '_H_' + H +
            '.dat', format='ascii')
        return pyUP_PHOT, pyUP_PM, AS_PHOT, AS_PM
    else:
        return pyUP_PHOT, pyUP_PM, UP_PHOT, UP_PM


def makePlot(
    N_UPMASK, tie_max, tie_min, m, win_PHOT, loss_PHOT, emp_PHOT,
        win_PM, loss_PM, emp_PM, flagasteca=False):
    """
    """
    category_names_PM = ['UPMASK_PM', 'Tie_PM', 'pyUPMASK_PM']
    category_names_PM_AS = ['ASteCA_PM', 'Tie_PM', 'pyUPMASK_PM']
    Results_PM = {
        'LSR': [loss_PM[0], emp_PM[0], win_PM[0]],
        'BSL': [loss_PM[1], emp_PM[1], win_PM[1]],
        'HMS': [loss_PM[2], emp_PM[2], win_PM[2]],
        r'MCC$_5$': [loss_PM[3], emp_PM[3], win_PM[3]],
        r'TPR$_5$': [loss_PM[4], emp_PM[4], win_PM[4]],
        r'PPV$_5$': [loss_PM[5], emp_PM[5], win_PM[5]],
        r'MCC$_9$': [loss_PM[6], emp_PM[6], win_PM[6]],
        r'TPR$_9$': [loss_PM[7], emp_PM[7], win_PM[7]],
        r'PPV$_9$': [loss_PM[8], emp_PM[8], win_PM[8]]
    }
    category_names_PHOT = ['UPMASK_PHOT', 'Tie_PHOT', 'pyUPMASK_PHOT']
    category_names_PHOT_AS = ['ASteCA_PHOT', 'Tie_PHOT', 'pyUPMASK_PHOT']
    Results_PHOT = {
        'LSR': [loss_PHOT[0], emp_PHOT[0], win_PHOT[0]],
        'BSL': [loss_PHOT[1], emp_PHOT[1], win_PHOT[1]],
        'HMS': [loss_PHOT[2], emp_PHOT[2], win_PHOT[2]],
        r'MCC$_5$': [loss_PHOT[3], emp_PHOT[3], win_PHOT[3]],
        r'TPR$_5$': [loss_PHOT[4], emp_PHOT[4], win_PHOT[4]],
        r'PPV$_5$': [loss_PHOT[5], emp_PHOT[5], win_PHOT[5]],
        r'MCC$_9$': [loss_PHOT[6], emp_PHOT[6], win_PHOT[6]],
        r'TPR$_9$': [loss_PHOT[7], emp_PHOT[7], win_PHOT[7]],
        r'PPV$_9$': [loss_PHOT[8], emp_PHOT[8], win_PHOT[8]]
    }

    # Combine photometric and PMS data
    category_names_comb = ['UPMASK', 'Tie', 'pyUPMASK']
    category_names_comb_as = ['ASteCA', 'Tie', 'pyUPMASK']
    Results_comb = {}
    for k, v in Results_PM.items():
        Results_comb[k] = np.array(v) + np.array(Results_PHOT[k])
    if flagasteca == True:
        return Results_PM, category_names_PM_AS
    else:
        return Results_PM, category_names_PM


def barsPlot(ax, results, category_names):
    """
    Parameters
    ----------
    results : dict
        A mapping from question labels to a list of answers per category.
        It is assumed all lists contain the same number of entries and that
        it matches the length of *category_names*.
    category_names : list of str
        The category labels.
    """
    labels = list(results.keys())
    data = np.array(list(results.values()))
    data_cum = data.cumsum(axis=1)
    category_colors = plt.get_cmap('RdYlBu')(
        np.linspace(0.15, 0.85, data.shape[1]))

    ax.invert_yaxis()
    ax.xaxis.set_visible(False)
    ax.set_xlim(0, np.sum(data, axis=1).max())

    for i, (colname, color) in enumerate(zip(category_names, category_colors)):
        widths = data[:, i]
        starts = data_cum[:, i] - widths
        ax.barh(labels, widths, left=starts, height=0.55,
                label=colname, color=color)
        xcenters = starts + widths / 2

        r, g, b, _ = color
        text_color = 'black'
        for y, (x, c) in enumerate(zip(xcenters, widths)):
            if c != 0:
                ax.text(x, y, str(int(c)), ha='center', va='center',
                        color=text_color, fontsize=12)
    ax.legend(ncol=len(category_names), bbox_to_anchor=(0, 1),
              loc='lower left', fontsize=11)


if __name__ == '__main__':
    main()