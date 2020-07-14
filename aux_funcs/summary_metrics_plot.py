
import numpy as np
import matplotlib.pyplot as plt
from metrics_plot import tie_min, tie_max, WinTieLoss, readTables


def main():
    """
    Make the summary plot that compares the different methods tested.
    """

    # Bad performers:
    # '75perc', 'marginNmemb_autoperc', 'marginC_2_autoperc',
    # 'GUMMprobs_autoperc', 'norm_GUMMprobs_autoperc', 'marginC_autoperc'
    # 'manualperc_1', 'optm_GUMM', 'autoperc', 'autoperc_5', 'autoperc_10',

    # This performs identical to 'autoperc_inner_GUMM3'
    # 'inner_GUMM_marginC'

    Hval = 'auto'  # 'symm', 'SR05')
    N_UPMASK = ("25", "50")

    # Folder where the files are located
    fold = "../TEST_SYNTH_CLUSTS/test_results/"

    for mode in ("VOR",): # VOR, GMM

        if mode == 'GMM':
            configs = (
                'autoperc_GMM', 'autoperc_GMM2', 'autoperc_GMM3',
                'autoperc_GMM4', 'minibatch_50', 'minibatch_vor')
        elif mode == "VOR":
            configs = (
                'autoperc_inner_GUMM', 'autoperc_inner_GUMM2',
                'autoperc_inner_GUMM3', 'autoperc_inner_GUMM4',
                'autoperc_inner_GUMM5', 'autoperc_inner_GUMM6',
                'autoperc_inner_GUMM7', 'voronoi_norm',
                'voronoi_dist', 'voronoi_2idx', 'voronoi_kmeans',
                'agglomerative_50', 'agglomerative_25')

        for NU in N_UPMASK:
            print(NU)

            winloss_rates = {}
            for m in configs:

                pyUP_PHOT, pyUP_PM, UP_PHOT, UP_PM = readTables(
                    fold, NU, Hval, m)

                CI_PM, CI_PHOT, win_PHOT, loss_PHOT, win_PM, loss_PM =\
                    WinTieLoss(
                        tie_max, tie_min, pyUP_PHOT, pyUP_PM, UP_PHOT,
                        UP_PM, 'summary')

                winloss_rates[m] = (
                    win_PHOT.sum() / loss_PHOT.sum(),
                    win_PM.sum() / loss_PM.sum(),
                    win_PHOT - loss_PHOT, win_PM - loss_PM)

            makePlot(fold, Hval, mode, NU, winloss_rates)


def makePlot(fold, H, mode, NU, winloss_rates):
    """
    Summary of the combined metrics for all the methods.
    """
    fig = plt.figure(figsize=(15, 10))

    min_y, max_y = np.inf, 0
    plt.subplot(211)
    for k, v in winloss_rates.items():
        win_phot, win_pm = (v[2] > 0.).sum(), (v[3] > 0.).sum()
        print("{} NU={}: PM {:.2f} ({} wins), PHOT {:.2f} ({} wins)".format(
            k, NU, v[1], win_pm, v[0], win_phot))

        plt.scatter(v[0], v[1], alpha=.7, s=50)
        plt.annotate(k, (v[0], v[1]))
        min_y = min(min_y, v[1] - .15)
        max_y = max(max_y, v[1] + .15)
    plt.xlabel("PHOT (W/L)")
    plt.ylabel("PM (W/L)")
    plt.ylim(min_y - .5, max_y + .5)

    labels = list(winloss_rates.keys())
    phot_x, pm_x = [], []
    for k, v in winloss_rates.items():
        win_phot, win_pm = (v[2] > 0.).sum(), (v[3] > 0.).sum()
        phot_x.append(win_phot)
        pm_x.append(win_pm)
    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    ax = plt.subplot(212)
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

    file_out = fold + 'plots/summary_{}_{}_{}.png'.format(H, mode, NU)
    fig.tight_layout()
    plt.savefig(file_out, dpi=150, bbox_inches='tight')


if __name__ == '__main__':
    main()
