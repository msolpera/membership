
import matplotlib.pyplot as plt
from metrics_plot import readTables


def main():
    """
    Plot the PPV vs the TPR (90%) emulating a Precision-Recall plot[1],
    although it is *not* a PR plot as no threshold is being varied.

    [1]: https://scikit-learn.org/stable/auto_examples/model_selection/
    plot_precision_recall.html
    """
    configs = ("Agglo", "kNNde", "Voron", 'MiniB', "KMean", "Gauss")
    Hval = 'auto'
    N_UPMASK = "25"

    for i, m in enumerate(configs):
        print(" ", m)
        pyUP_PHOT, pyUP_PM, UP_PHOT, UP_PM = readTables(N_UPMASK, Hval, m)

        plt.subplot(int("23" + str(i + 1)))
        plt.title(m)
        plt.scatter(
            pyUP_PM['TPR_9'], pyUP_PM['PPV_9'], marker="o", c=pyUP_PM['CI'],
            label='PM', alpha=.5)
        plt.scatter(
            pyUP_PHOT['TPR_9'], pyUP_PHOT['PPV_9'], marker="s",
            c=pyUP_PHOT['CI'], label='PHOT', alpha=.5)
        if i == 0:
            plt.legend()
            plt.colorbar()
        plt.xlabel("Recall (TPR_9)")
        plt.ylabel("Precision (PPV_9)")
        plt.xlim(0., 1.05)
        plt.ylim(0., 1.05)
    plt.show()

    # Now UPMASK results
    plt.title("UPMASK")
    plt.scatter(
        UP_PM['TPR_9'], UP_PM['PPV_9'], marker="o", c=UP_PM['CI'],
        label='PM', alpha=.5)
    plt.scatter(
        UP_PHOT['TPR_9'], UP_PHOT['PPV_9'], marker="s", c=UP_PHOT['CI'],
        label='PHOT', alpha=.5)
    plt.legend()
    plt.colorbar()
    plt.xlabel("Recall (TPR_9)")
    plt.ylabel("Precision (PPV_9)")
    plt.xlim(0., 1.05)
    plt.ylim(0., 1.05)
    plt.show()


if __name__ == '__main__':
    main()
