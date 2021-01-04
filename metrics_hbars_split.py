
import numpy as np
import matplotlib.pyplot as plt
from auxFuncs import tie_max, tie_min, readTables, WinTieLoss


def main():
    """
    For each clustering method plot the horizontal bar plots for all the
    metrics, separated into PM, PHOT, and combined results.
    """
    configs = ("Agglo", "Gauss", "KMean", "kNNde", 'MiniB', "Voron")
    Hval = 'auto'  # 'symm', 'SR05'

    # Defines which UPMASK results to use
    # Cantat-Gaudin results
    # UP_alg, N_UPMASK = "CG_", "15" # "25"
    # Original UPMASK results
    UP_alg, N_UPMASK = "", "25"

    for m in configs:
        print(m)
        pyUP_PHOT, pyUP_PM, UP_PHOT, UP_PM = readTables(
            N_UPMASK, Hval, m, UP_alg)
        win_PHOT, loss_PHOT, emp_PHOT, win_PM, loss_PM, emp_PM =\
            WinTieLoss(pyUP_PHOT, pyUP_PM, UP_PHOT, UP_PM, 'metrics_bars')
        makePlot(
            N_UPMASK, tie_max, tie_min, m, win_PHOT, loss_PHOT,
            emp_PHOT, win_PM, loss_PM, emp_PM)


def makePlot(
    N_UPMASK, tie_max, tie_min, m, win_PHOT, loss_PHOT, emp_PHOT,
        win_PM, loss_PM, emp_PM):
    """
    """
    category_names_PM = ['Loss_PM', 'Tie_PM', 'Win_PM']
    Results_PM = {
        'LSR': [loss_PM[0], emp_PM[0], win_PM[0]],
        'BSL': [loss_PM[1], emp_PM[1], win_PM[1]],
        'HMS': [loss_PM[2], emp_PM[2], win_PM[2]],
        'MCC_5': [loss_PM[3], emp_PM[3], win_PM[3]],
        'TPR_5': [loss_PM[4], emp_PM[4], win_PM[4]],
        'PPV_5': [loss_PM[5], emp_PM[5], win_PM[5]],
        'MCC_9': [loss_PM[6], emp_PM[6], win_PM[6]],
        'TPR_9': [loss_PM[7], emp_PM[7], win_PM[7]],
        'PPV_9': [loss_PM[8], emp_PM[8], win_PM[8]]
    }
    category_names_PHOT = ['Loss_PHOT', 'Tie_PHOT', 'Win_PHOT']
    Results_PHOT = {
        'LSR': [loss_PHOT[0], emp_PHOT[0], win_PHOT[0]],
        'BSL': [loss_PHOT[1], emp_PHOT[1], win_PHOT[1]],
        'HMS': [loss_PHOT[2], emp_PHOT[2], win_PHOT[2]],
        'MCC_5': [loss_PHOT[3], emp_PHOT[3], win_PHOT[3]],
        'TPR_5': [loss_PHOT[4], emp_PHOT[4], win_PHOT[4]],
        'PPV_5': [loss_PHOT[5], emp_PHOT[5], win_PHOT[5]],
        'MCC_9': [loss_PHOT[6], emp_PHOT[6], win_PHOT[6]],
        'TPR_9': [loss_PHOT[7], emp_PHOT[7], win_PHOT[7]],
        'PPV_9': [loss_PHOT[8], emp_PHOT[8], win_PHOT[8]]
    }

    # Combine photometric and PMS data
    category_names_comb = ['Loss', 'Tie', 'Win']
    Results_comb = {}
    for k, v in Results_PM.items():
        Results_comb[k] = np.array(v) + np.array(Results_PHOT[k])

    plt.figure(figsize=(10, 12))
    plt.suptitle("Tie range: [{}, {}]".format(tie_min, tie_max), y=.95)
    ax = plt.subplot(311)
    barsPlot(ax, Results_PM, category_names_PM)
    ax = plt.subplot(312)
    barsPlot(ax, Results_PHOT, category_names_PHOT)
    ax = plt.subplot(313)
    barsPlot(ax, Results_comb, category_names_comb)
    file_out = 'plots/' + '{}_{}.png'.format(m, N_UPMASK)
    # plt.savefig(file_out, dpi=150, bbox_inches='tight')
    plt.show()


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
    category_colors = plt.get_cmap('RdYlGn')(
        np.linspace(0.15, 0.85, data.shape[1]))

    ax.invert_yaxis()
    ax.xaxis.set_visible(False)
    ax.set_xlim(0, np.sum(data, axis=1).max())

    for i, (colname, color) in enumerate(zip(category_names, category_colors)):
        widths = data[:, i]
        starts = data_cum[:, i] - widths
        ax.barh(labels, widths, left=starts, height=0.5,
                label=colname, color=color)
        xcenters = starts + widths / 2

        r, g, b, _ = color
        # text_color = 'white' if r * g * b < 0.5 else 'darkgrey'
        text_color = 'black'
        for y, (x, c) in enumerate(zip(xcenters, widths)):
            ax.text(x, y, str(int(c)), ha='center', va='center',
                    color=text_color, fontsize=12)
    ax.legend(ncol=len(category_names), bbox_to_anchor=(0, 1),
              loc='lower left', fontsize=11)


# def readFiles(ruta=Path.cwd()):
#     """
#     Read files from the input folder
#     """
#     return [arch.name for arch in Path('input').iterdir() if arch.is_file()]


if __name__ == '__main__':
    main()
