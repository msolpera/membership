
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from auxFuncs import tie_min, tie_max, WinTieLoss, readTables


def main():
    """
    Summary plot: win vs loss for PM & PHOT, and bar plots for each metric
    showing how many were won/lost
    """
    # Defines which results to use
    # # Cantat-Gaudin results
    UP_alg, N_UPMASK = "CG_", "15"  # "25"
    # Original UPMASK results
    # UP_alg, N_UPMASK = "", "25"

    Hval = 'auto'  # 'symm', 'SR05')
    configs = ("Agglo", "Gauss", "KMean", "kNNde", 'MiniB', "Voron")

    winloss_rates = {}
    for m in configs:

        pyUP_PHOT, pyUP_PM, UP_PHOT, UP_PM = readTables(
            N_UPMASK, Hval, m, UP_alg)

        CI_PM, CI_PHOT, win_PHOT, loss_PHOT, win_PM, loss_PM =\
            WinTieLoss(pyUP_PHOT, pyUP_PM, UP_PHOT, UP_PM, 'summary')

        delta_PHOT, delta_PM = win_PHOT - loss_PHOT, win_PM - loss_PM
        delta_sum_PHOT = win_PHOT.sum() - loss_PHOT.sum()
        delta_sum_PM = win_PM.sum() - loss_PM.sum()

        N_tot_PM = CI_PM.size * win_PM.size
        N_tot_PHOT = CI_PHOT.size * win_PHOT.size

        winloss_rates[m] = (
            100 * (loss_PHOT.sum() / N_tot_PHOT),
            100 * (loss_PM.sum() / N_tot_PM),
            100 * (win_PHOT.sum() / N_tot_PHOT),
            100 * (win_PM.sum() / N_tot_PM),
            100 * (delta_sum_PHOT / N_tot_PHOT),
            100 * (delta_sum_PM / N_tot_PM),
            delta_PHOT, delta_PM)

    makePlot(Hval, N_UPMASK, tie_max, tie_min, winloss_rates)


def makePlot(H, NU, tie_max, tie_min, winloss_rates):
    """
    Summary of the combined metrics for all the methods.
    """
    fig = plt.figure(figsize=(15, 10))

    G = gridspec.GridSpec(3, 2)
    plt.suptitle("Tie range: [{}, {}]".format(tie_min, tie_max), y=.95)

    min_y, max_y, min_x, max_x = np.inf, 0, np.inf, 0.
    plt.subplot(G[0, 0])
    plt.title("Percentage lost to UPMASK")
    for k, v in winloss_rates.items():
        win_phot, win_pm = (v[-2] > 0.).sum(), (v[-1] > 0.).sum()
        print("\n" + k, "NU=", NU)
        print("PM   W {:.2f}, T {:.2f}, L {:.2f} ({} wins)".format(
            v[3], 100. - (v[3] + v[1]), v[1], win_pm))
        print("PHOT W {:.2f}, T {:.2f}, L {:.2f} ({} wins)".format(
            v[2], 100. - (v[2] + v[0]), v[0], win_phot))
        mrk = 'o'
        if 'GMM' in k or 'minibatch' in k:
            mrk = 's'
        plt.scatter(v[0], v[1], marker=mrk, alpha=.5, s=150)
        plt.annotate(k, (v[0], v[1]))
        min_y = min(min_y, v[1] - .15)
        max_y = max(max_y, v[1] + .15)
        min_x = min(min_x, v[0] - .15)
        max_x = max(max_x, v[0] + .15)
    plt.xlabel("PHOT L%")
    plt.ylabel("PM L%")
    #
    plt.subplot(G[0, 1])
    plt.title("Percentage won to UPMASK")
    for k, v in winloss_rates.items():
        mrk = 'o'
        if 'GMM' in k or 'minibatch' in k:
            mrk = 's'
        plt.scatter(v[2], v[3], marker=mrk, alpha=.5, s=150)
        plt.annotate(k, (v[2], v[3]))
        min_y = min(min_y, v[3] - .15)
        max_y = max(max_y, v[3] + .15)
        min_x = min(min_x, v[2] - .15)
        max_x = max(max_x, v[2] + .15)
    plt.xlabel("PHOT W%")
    plt.ylabel("PM W%")

    min_y, max_y, min_x, max_x = np.inf, 0, np.inf, 0.
    plt.subplot(G[1, 0:2])
    for k, v in winloss_rates.items():
        win_phot, win_pm = (v[-2] > 0.).sum(), (v[-1] > 0.).sum()
        print(
            "{} NU={}: PM {:.2f} ({} wins), PHOT {:.2f} ({} wins)".format(
                k, NU, v[5], win_pm, v[4], win_phot))
        mrk = 'o'
        if 'GMM' in k or 'minibatch' in k:
            mrk = 's'
        plt.scatter(v[4], v[5], marker=mrk, alpha=.5, s=150)
        plt.annotate(k, (v[4], v[5]))
        min_y = min(min_y, v[5] - .15)
        max_y = max(max_y, v[5] + .15)
        min_x = min(min_x, v[4] - .15)
        max_x = max(max_x, v[4] + .15)
    plt.axvline(0., ls=':')
    plt.axhline(0., ls=':')
    # plt.axvspan(-10, 0, facecolor='r', alpha=0.25)
    plt.xlabel("PHOT (W-L)%")
    plt.ylabel("PM (W-L)%")
    plt.xlim(min_x - 5, max_x + 5)
    plt.ylim(min_y - 5, max_y + 5)

    ax = plt.subplot(G[2, 0:2])
    labels = list(winloss_rates.keys())
    phot_x, pm_x = [], []
    for k, v in winloss_rates.items():
        win_phot, win_pm = (v[-2] > 0.).sum(), (v[-1] > 0.).sum()
        phot_x.append(win_phot)
        pm_x.append(win_pm)
    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    h_phot = ax.bar(x - width / 2, phot_x, width, label='PHOT')
    h_pm = ax.bar(x + width / 2, pm_x, width, label='PM')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)

    def autolabel(rects):
        for rect in rects:
            height = rect.get_height()
            ax.annotate('{}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')

    autolabel(h_phot)
    autolabel(h_pm)
    plt.xticks(rotation=45)
    plt.ylim(0, 10)
    plt.legend()

    plt.show()

    # file_out = 'plots/summary_{}_{}.png'.format(H, NU)
    # fig.tight_layout()
    # plt.savefig(file_out, dpi=150, bbox_inches='tight')


if __name__ == '__main__':
    main()