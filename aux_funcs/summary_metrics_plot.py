
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
    # 'manualperc_1'

    # This performs identical to 'autoperc_inner_GUMM3'
    # 'inner_GUMM_marginC'

    mode = (
        'autoperc', 'autoperc_5', 'autoperc_10', 'autoperc_inner_GUMM',
        'autoperc_inner_GUMM2', 'autoperc_inner_GUMM3',
        'autoperc_inner_GUMM4', 'autoperc_inner_GUMM5', 'autoperc_GMM',
        'autoperc_GMM2', 'autoperc_GMM3', 'optm_GUMM')
    Hval = ('auto',)  # 'symm', 'SR05')

    # Folder where the files are located
    fold = "../TEST_SYNTH_CLUSTS/test_results/"

    for H in Hval:
        print(H)
        winloss_rates = {}
        for m in mode:
            print(" ", m)
            pyUP_PHOT, pyUP_PM, UP_PHOT, UP_PM = readTables(fold, H, m)

            CI_PM, CI_PHOT, win_PHOT, loss_PHOT, win_PM, loss_PM = WinTieLoss(
                tie_max, tie_min, pyUP_PHOT, pyUP_PM, UP_PHOT, UP_PM,
                'summary')

            winloss_rates[m] = (
                win_PHOT.sum() / loss_PHOT.sum(), win_PM.sum() / loss_PM.sum(),
                win_PHOT - loss_PHOT, win_PM - loss_PM)

        makePlot(fold, H, winloss_rates)


def makePlot(fold, H, winloss_rates):
    """
    Summary of the combined metrics for all the methods.
    """
    fig = plt.figure(figsize=(15, 10))

    min_y, max_y = np.inf, 0
    plt.subplot(211)
    i = 0
    for k, v in winloss_rates.items():
        yoff = -.1 if (i % 2) == 0 else .07
        i += 1
        xr, yr = np.random.uniform(.015, .02, 2)
        plt.scatter(v[0] + xr, v[1] + yr, alpha=.7, s=50)
        # plt.annotate(k, (v[0] - .025, v[1] + yoff))
        plt.annotate(k, (v[0], v[1] + yoff))
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

    file_out = fold + 'plots/summary_{}.png'.format(H)
    fig.tight_layout()
    plt.savefig(file_out, dpi=150, bbox_inches='tight')


if __name__ == '__main__':
    main()
