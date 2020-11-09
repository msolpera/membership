
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from metrics_vert_bars import tie_min, tie_max, WinTieLoss, readTables


def main():
    """
    """
    Hval = 'auto'  # 'symm', 'SR05')
    N_UPMASK = "25"  # "50")
    configs = ("Voron", "kNNde", "Agglo", 'MiniB', "KMean", "Gauss")
    col_labels = ('VOR', 'KNN', 'AGG', 'MBK', 'KMS', 'GMM')
    metrics = ["LSR", "BSL", "HMS", r"MCC$_5$", r"TPR$_5$", r"PPV$_5$",
               r"MCC$_9$", r"TPR$_9$", r"PPV$_9$"]

    winloss_rates = [[], [], []]
    for m in configs:

        pyUP_PHOT, pyUP_PM, UP_PHOT, UP_PM = readTables(N_UPMASK, Hval, m)

        CI_PM, CI_PHOT, win_PHOT, loss_PHOT, win_PM, loss_PM =\
            WinTieLoss(
                tie_max, tie_min, pyUP_PHOT, pyUP_PM, UP_PHOT,
                UP_PM, 'summary')

        # PM
        delta = win_PM - loss_PM
        N_tot = CI_PM.size
        winloss_rates[0].append(100 * (delta / N_tot))

        # PHOT
        delta = win_PHOT - loss_PHOT
        N_tot = CI_PHOT.size
        winloss_rates[1].append(100 * (delta / N_tot))

        # Delta for pyUPMASK
        delta = (win_PHOT + win_PM) - (loss_PHOT + loss_PM)
        # # Delta for UPMASK
        # delta = (loss_PHOT + loss_PM) - (win_PHOT + win_PM)

        N_tot = CI_PM.size + CI_PHOT.size
        winloss_rates[2].append(100 * (delta / N_tot))

    titl = ["PM", "PHOT", "Combined"]
    for i, win_loss in enumerate(winloss_rates):
        print(titl[i])
        fig, ax = plt.subplots(figsize=(7, 7))
        # plt.title(titl[i])
        # (metrics, methods)
        matrx_vals = np.array(win_loss).T
        im, cbar = heatmap(matrx_vals, metrics, col_labels, ax=ax,
                           cmap="RdYlGn", cbarlabel="(W-L)%") # YlGn
        annotate_heatmap(im, valfmt="{x:.0f}")
        # plt.show()
        file_out = 'plots/matrix_{}.png'.format(titl[i])
        fig.tight_layout()
        plt.savefig(file_out, dpi=300, bbox_inches='tight')


def heatmap(
    data, row_labels, col_labels, ax=None, cbar_kw={}, cbarlabel="",
        **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, vmin=-100., vmax=100, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1] + 1) - .5, minor=True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - .5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(
    im, data=None, valfmt="{x:.2f}", textcolors=("black", "black"),
        threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max()) / 2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = round(data[i, j], 0)
            # This line prevents '-0' values
            val = 0. if val == 0. else val
            if val < 0.:
                clr = 'black' #'red'
            else:
                clr = textcolors[int(im.norm(val) > threshold)]
            kw.update(color=clr)
            text = im.axes.text(j, i, valfmt(val, None), **kw)
            texts.append(text)

    return texts


if __name__ == '__main__':
    main()
