
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from metrics_vert_bars import tie_min, tie_max, WinTieLoss, readTables


def main():
    """
    Make the summary plot that compares the different methods tested.
    """
    Hval = 'auto'  # 'symm', 'SR05')
    N_UPMASK = "25"  # "50")
    configs = ("Voron", "kNNde", "Agglo", 'MiniB', "KMean", "Gauss")

    winloss_rates = {}
    for m in configs:

        pyUP_PHOT, pyUP_PM, UP_PHOT, UP_PM = readTables(
            N_UPMASK, Hval, m)

        CI_PM, CI_PHOT, win_PHOT, loss_PHOT, win_PM, loss_PM =\
            WinTieLoss(
                tie_max, tie_min, pyUP_PHOT, pyUP_PM, UP_PHOT,
                UP_PM, 'summary')

        win_sum = win_PM.sum() + win_PHOT.sum()
        loss_sum = loss_PM.sum() + loss_PHOT.sum()

        N_tot_PM = CI_PM.size * win_PM.size
        N_tot_PHOT = CI_PHOT.size * win_PHOT.size
        N_tot = N_tot_PM + N_tot_PHOT

        winloss_rates[m] = (
            100 * (win_sum / N_tot),
            100 * (loss_sum / N_tot),
            100 * (N_tot - win_sum - loss_sum) / N_tot)

    makePlot(Hval, N_UPMASK, tie_max, tie_min, winloss_rates)


def makePlot(H, NU, tie_max, tie_min, winloss_rates):
    """
    Summary of the combined metrics for all the methods.
    """
    fig = plt.figure(figsize=(10, 10))

    # G = gridspec.GridSpec(2, 2)
    plt.suptitle("pyUPMASK vs UPMASK", y=.9, fontsize=10)
    plt.grid(ls=':', c='grey', lw=.5)

    min_y, max_y, min_x, max_x = np.inf, 0, np.inf, 0.
    for k, v in winloss_rates.items():
        plt.scatter(v[0], v[1], marker='o', s=250, alpha=.7)
        names = {"Agglo": 'AGG', "Gauss": 'GMM', "KMean": 'KMS',
                 "kNNde": 'KNN', 'MiniB': 'MBK', "Voron": 'VOR'}
        print(names[k], v[0], v[1])
        plt.annotate(names[k] + " ({:.0f}%)".format(
            v[2]), (v[0] + .5, v[1] + .5), fontsize=10)
        min_y = min(min_y, v[1] - .15)
        max_y = max(max_y, v[1] + .15)
        min_x = min(min_x, v[0] - .15)
        max_x = max(max_x, v[0] + .15)

    plt.xlabel("WIN %", fontsize=10)
    plt.ylabel("LOSS %", fontsize=10)
    plt.xlim(min_x - 5, max_x + 5)
    plt.ylim(min_y - 5, max_y + 5)
    # plt.xlim(50, 100)
    # plt.ylim(0, 50)

    plt.show()

    # file_out = 'plots/summary_{}_{}.png'.format(H, NU)
    # fig.tight_layout()
    # plt.savefig(file_out, dpi=150, bbox_inches='tight')


if __name__ == '__main__':
    main()
