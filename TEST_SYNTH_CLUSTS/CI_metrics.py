
import matplotlib.pyplot as plt
from metrics_plot import tie_min, tie_max, WinTieLoss, readTables


def main(plot_raw_metrics=False):
    """
    Plot the metrics versus CI for all the clusters.
    """
    configs = ("Voron", "kNNde", "Agglo", 'MiniB', "KMean", "Gauss")
    Hval = 'auto'  # 'symm', 'SR05'
    N_UPMASK = "25"  # "50"

    for m in configs:
        print(" ", m)
        pyUP_PHOT, pyUP_PM, UP_PHOT, UP_PM = readTables(N_UPMASK, Hval, m)

        CI_PM, CI_PHOT, comp_PHOT, comp_PM, pyU_PHOT, pyU_PM = WinTieLoss(
            tie_max, tie_min, pyUP_PHOT, pyUP_PM, UP_PHOT, UP_PM, 'CI')

        makePlot(
            tie_max, tie_min, m, CI_PM, CI_PHOT, comp_PHOT,
            comp_PM, pyU_PHOT, pyU_PM, plot_raw_metrics)


def makePlot(
    tie_max, tie_min, m, CI_PM, CI_PHOT, comp_PHOT, comp_PM,
        pyU_PHOT, pyU_PM, plot_raw_metrics):
    """
    """
    metrics = (
        'LSR', 'BSL', 'HMS', 'MCC_5', 'TPR_5', 'PPV_5', 'MCC_9', 'TPR_9',
        'PPV_9')

    fig = plt.figure(figsize=(20, 10))
    CIPlot(metrics, tie_min, tie_max, CI_PM, comp_PM, CI_PHOT, comp_PHOT)
    file_out = 'plots/' + 'CI_delta_{}.png'.format(m)
    fig.tight_layout()
    plt.savefig(file_out, dpi=150, bbox_inches='tight')
    # plt.show()

    if plot_raw_metrics:
        fig = plt.figure(figsize=(20, 10))
        rawMetricsPlot(metrics, CI_PM, CI_PHOT, pyU_PHOT, pyU_PM)
        file_out = 'plots/' + 'CI_{}.png'.format(m)
        fig.tight_layout()
        plt.savefig(file_out, dpi=150, bbox_inches='tight')
        # plt.show()


def CIPlot(metrics, tie_min, tie_max, CI_PM, comp_PM, CI_PHOT, comp_PHOT):
    """
    """
    for i, met in enumerate(metrics):
        ax = plt.subplot(int("33" + str(i + 1)))
        plt.title(met)
        plt.axhline(0., c='k', ls=':')
        plt.scatter(
            CI_PM, comp_PM[i], c='k', alpha=.5, label="PM", zorder=3)
        plt.scatter(
            CI_PHOT, comp_PHOT[i], c='b', marker='^', alpha=.5, label="PHOT",
            zorder=3)
        ymin, ymax = ax.get_ylim()

        ax.axhspan(tie_min, tie_max, alpha=.1, color='yellow')
        ax.axhspan(-10, tie_min, alpha=.1, color='red')
        ax.axhspan(tie_max, 10, alpha=.1, color='green')
        plt.ylim(ymin, ymax)

        if met in ('LSR', 'MCC_5', 'MCC_9'):
            plt.ylabel("pyUPMASK - UPMASK")
        if met in ('MCC_9', 'TPR_9', 'PPV_9'):
            plt.xlabel("CI (log)")
        ax.set_xscale('log')
        if i == 0:
            plt.legend()


def rawMetricsPlot(metrics, CI_PM, CI_PHOT, pyU_PHOT, pyU_PM):
    """
    """
    for i, met in enumerate(metrics):
        ax = plt.subplot(int("33" + str(i + 1)))
        plt.title(met)
        plt.axhline(0., c='k', ls=':')
        plt.scatter(
            CI_PM, pyU_PM[i], c='k', alpha=.5, label="PM", zorder=3)
        plt.scatter(
            CI_PHOT, pyU_PHOT[i], c='b', marker='^', alpha=.5, label="PHOT",
            zorder=3)

        if met in ('LSR', 'MCC_5', 'MCC_9'):
            plt.ylabel(r"$metric$")
        if met in ('MCC_9', 'TPR_9', 'PPV_9'):
            plt.xlabel("CI (log)")
        ax.set_xscale('log')
        if i == 0:
            plt.legend()


if __name__ == '__main__':
    main()
