
import numpy as np
import matplotlib.pyplot as plt
from metrics_vert_bars import tie_min, tie_max, WinTieLoss, readTables,\
    barsPlot


def main():
    """
    For all the clustering methods plot the horizontal bar plots for all the
    metrics, separated into PM, PHOT, and combined results.
    """
    configs = ("Voron", "kNNde", "Agglo", "MiniB", "KMean", "Gauss")
    Hval = 'auto'  # 'symm', 'SR05'
    N_UPMASK = "25"  # "50"

    results = []
    for i, m in enumerate(configs):
        print(m)
        pyUP_PHOT, pyUP_PM, UP_PHOT, UP_PM = readTables(
            N_UPMASK, Hval, m)
        win_PHOT, loss_PHOT, emp_PHOT, win_PM, loss_PM, emp_PM =\
            WinTieLoss(tie_max, tie_min, pyUP_PHOT, pyUP_PM, UP_PHOT,
                       UP_PM, 'metrics_bars')
        results.append(combResults(
            N_UPMASK, tie_max, tie_min, m, win_PHOT, loss_PHOT,
            emp_PHOT, win_PM, loss_PM, emp_PM))

    # Make final plot
    fig = plt.figure(figsize=(13, 10))
    conf = ("VOR", "KNN", "AGG", "MBK", "KMS", "GMM")
    for j, l in enumerate(conf):
        ax = plt.subplot(3, 2, j + 1)
        barsPlot(ax, results[j][0], results[j][1])
        plt.title(l, fontsize=10)
        ax.invert_xaxis()
        if j != 0:
            ax.legend().remove()
    # fig.tight_layout()
    # plt.savefig('metrics_bars.png', dpi=300, bbox_inches='tight')
    plt.show()


def combResults(
    N_UPMASK, tie_max, tie_min, m, win_PHOT, loss_PHOT, emp_PHOT,
        win_PM, loss_PM, emp_PM):
    """
    """
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
    category_names_comb = ['Loss', 'Tie', 'Win']
    Results_comb = {}
    for k, v in Results_PM.items():
        Results_comb[k] = np.array(v) + np.array(Results_PHOT[k])

    return Results_comb, category_names_comb


if __name__ == '__main__':
    main()
