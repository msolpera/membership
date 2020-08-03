
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
    # 'voronoi_norm', 'autoperc_inner_GUMM', 'autoperc_inner_GUMM2'
    # 'autoperc_inner_GUMM5', 'voronoi_kmeans', 'autoperc_inner_GUMM7'
    # , 'voronoi_2idx'  'autoperc_inner_GUMM6', 'autoperc_inner_GUMM3'
    #  'voronoi_dist', 'agglomerative_KDE', 'kNN_10_10'
    # 'voronoi_kde_p_25', , 'minibatch_vor'

    # This performs identical to 'autoperc_inner_GUMM3'
    # 'inner_GUMM_marginC'
    # 'agglomerative_KDEpy' similar to 'agglomerative_KDE'
    # 'autoperc_inner_GMM4' similar to 'autoperc_inner_GMM3'
    # 'autoperc_GMM' similar to 'autoperc_GMM2'
    # 'kNN_25_mean_25' similar to 'kNN_25_25'
    # 'kNN_10_25', 'kNN_50_25', 'kNN_25_auto', 'kNN_10_50', 'kNN_25_50'
    # superseded by 'voronoi_kde_p_25': 'autoperc_inner_GUMM2',
    # 'autoperc_inner_GUMM6', 'voronoi_kde_p_50'
    # Good PM, very bad PHOT: 'agglomerative_KDEpy'
    # Superseded by 'kNN_50_50_kde': 'kNN_50_50'
    # Superseded by 'kNN_250_25_kde': 'kNN_25_25'
    # Superseded by 'agglomerative_25_kde_p': 'agglomerative_25',
    # 'agglomerative_50_kde_p'
    # Superseded by 'voronoi_newcents_25': 'voronoi_newcents_50',
    # 'voronoi_newcents_15'
    # Superseded by 'kNN_25_25_kde': 'kNN_50_50_kde'

    Hval = 'auto'  # 'symm', 'SR05')
    N_UPMASK = "25"  # "50")

    # Folder where the files are located
    fold = "../TEST_SYNTH_CLUSTS/test_results/"

    configs = (
        'autoperc_GMM2', 'autoperc_GMM3', 'minibatch_50',
        'kNN_25_25_kde', 'agglomerative_25_kde_p',
        'voronoi_newcents_25', 'knn_test')

    # for NU in N_UPMASK:
    #     print(NU)

    winloss_rates = {}
    for m in configs:

        pyUP_PHOT, pyUP_PM, UP_PHOT, UP_PM = readTables(
            fold, N_UPMASK, Hval, m)

        CI_PM, CI_PHOT, win_PHOT, loss_PHOT, win_PM, loss_PM =\
            WinTieLoss(
                tie_max, tie_min, pyUP_PHOT, pyUP_PM, UP_PHOT,
                UP_PM, 'summary')

        delta_PHOT, delta_PM = win_PHOT - loss_PHOT, win_PM - loss_PM
        delta_sum_PHOT = win_PHOT.sum() - loss_PHOT.sum()
        delta_sum_PM = win_PM.sum() - loss_PM.sum()
        N_tot_PM = CI_PM.size * win_PM.size
        N_tot_PHOT = CI_PHOT.size * win_PHOT.size
        winloss_rates[m] = (
            100 * (delta_sum_PHOT / N_tot_PHOT),
            100 * (delta_sum_PM / N_tot_PM),
            delta_PHOT, delta_PM)

    makePlot(fold, Hval, N_UPMASK, tie_max, tie_min, winloss_rates)


def makePlot(fold, H, NU, tie_max, tie_min, winloss_rates):
    """
    Summary of the combined metrics for all the methods.
    """
    fig = plt.figure(figsize=(15, 10))
    plt.suptitle("Tie range: [{}, {}]".format(tie_min, tie_max), y=.95)

    min_y, max_y, min_x, max_x = np.inf, 0, np.inf, 0.
    plt.subplot(211)
    for k, v in winloss_rates.items():
        win_phot, win_pm = (v[2] > 0.).sum(), (v[3] > 0.).sum()
        print("{} NU={}: PM {:.2f} ({} wins), PHOT {:.2f} ({} wins)".format(
            k, NU, v[1], win_pm, v[0], win_phot))

        mrk = 'o'
        if 'GMM' in k or 'minibatch' in k:
            mrk = 's'
        plt.scatter(v[0], v[1], marker=mrk, alpha=.5, s=150)
        plt.annotate(k, (v[0], v[1]))
        min_y = min(min_y, v[1] - .15)
        max_y = max(max_y, v[1] + .15)
        min_x = min(min_x, v[0] - .15)
        max_x = max(max_x, v[0] + .15)

    plt.axvline(0., ls=':')
    plt.axhline(0., ls=':')
    # plt.axvspan(-10, 0, facecolor='r', alpha=0.25)

    plt.xlabel("PHOT (W-L)%")
    plt.ylabel("PM (W-L)%")
    plt.xlim(min_x - 5, max_x + 5)
    # plt.ylim(min(0., min_y - 5), max_y + 5)
    plt.ylim(min_y - 5, max_y + 5)

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

    plt.show()

    # file_out = fold + 'plots/summary_{}_{}.png'.format(H, NU)
    # fig.tight_layout()
    # plt.savefig(file_out, dpi=150, bbox_inches='tight')


if __name__ == '__main__':
    main()
