from pathlib import Path
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
# from matplotlib.pyplot import figure
# from astropy.io import ascii


def main(summaryPlotFlag=True):
    """
    Plot the statistics obtained after applying the getMetrics() script to
    the pyUPMASK and UPMASK results.
    """

    # Define the "tie" range.
    tie_max = .05
    tie_min = -1. * tie_max

    # Bad performers:
    # '75perc', 'marginNmemb_autoperc', 'marginC_2_autoperc',
    # 'GUMMprobs_autoperc'
    mode = (
        'autoperc', 'marginC_autoperc', 'autoperc_5', 'autoperc_10',
        'norm_GUMMprobs_autoperc', 'manualperc_1', 'autoperc_inner_GUMM',
        'autoperc_inner_GUMM2', 'autoperc_inner_GUMM3', 'inner_GUMM_marginC',
        'autoperc_inner_GUMM4', 'autoperc_inner_GUMM5')
    # mode = ('inner_GUMM_marginC',)
    Hval = ('auto', 'symm', 'SR05')

    # Folder where the files are located
    fold = "../TEST_SYNTH_CLUSTS/test_results/"

    for H in Hval:
        print(H)
        winloss_rates = {}
        for m in mode:
            print(" ", m)
            pyUP_PHOT = Table.read(
                fold + 'metrics_PHOT_{}_H_{}.dat'.format(m, H), format='ascii')
            pyUP_PM = Table.read(
                fold + 'metrics_PM_{}_H_{}.dat'.format(m, H), format='ascii')
            UP_PHOT = Table.read(
                fold + 'metrics_UP-PHOT_H_' + H + '.dat', format='ascii')
            UP_PM = Table.read(
                fold + 'metrics_UP-PM_H_' + H + '.dat', format='ascii')

            CI_PM, CI_PHOT, comp_PHOT, comp_PM, win_PHOT, loss_PHOT, emp_PHOT,\
                win_PM, loss_PM, emp_PM = WinTieLoss(
                    tie_max, tie_min, pyUP_PHOT, pyUP_PM, UP_PHOT, UP_PM)

            # Used by summaryPlot()
            winloss_rates[m] = (
                win_PHOT.sum() / loss_PHOT.sum(), win_PM.sum() / loss_PM.sum(),
                win_PHOT - loss_PHOT, win_PM, loss_PM)

            if summaryPlotFlag is False:
                makePlot(
                    fold, tie_max, tie_min, H, m, CI_PM, CI_PHOT, comp_PHOT,
                    comp_PM, win_PHOT, loss_PHOT, emp_PHOT, win_PM, loss_PM,
                    emp_PM)

        if summaryPlotFlag:
            summaryPlot(fold, H, winloss_rates)


def WinTieLoss(tie_max, tie_min, pyUP_PHOT, pyUP_PM, UP_PHOT, UP_PM):
    """
    """
    CI_PM, CI_PHOT = pyUP_PM['CI'], pyUP_PHOT['CI']

    LSR_pUPH, BSL_pUPH, HMS_pUPH, MCC_5_pUPH, TPR_5_pUPH, PPV_5_pUPH,\
        MCC_9_pUPH, TPR_9_pUPH, PPV_9_pUPH = pyUP_PHOT['LSR'],\
        pyUP_PHOT['BSL'], pyUP_PHOT['HMS'], pyUP_PHOT['MCC_5'],\
        pyUP_PHOT['TPR_5'], pyUP_PHOT['PPV_5'], pyUP_PHOT['MCC_9'],\
        pyUP_PHOT['TPR_9'], pyUP_PHOT['PPV_9']

    LSR_pUPM, BSL_pUPM, HMS_pUPM, MCC_5_pUPM, TPR_5_pUPM, PPV_5_pUPM,\
        MCC_9_pUPM, TPR_9_pUPM, PPV_9_pUPM = pyUP_PM['LSR'], pyUP_PM['BSL'],\
        pyUP_PM['HMS'], pyUP_PM['MCC_5'], pyUP_PM['TPR_5'], pyUP_PM['PPV_5'],\
        pyUP_PM['MCC_9'], pyUP_PM['TPR_9'], pyUP_PM['PPV_9']

    LSR_UPH, BSL_UPH, HMS_UPH, MCC_5_UPH, TPR_5_UPH, PPV_5_UPH, MCC_9_UPH,\
        TPR_9_UPH, PPV_9_UPH = UP_PHOT['LSR'], UP_PHOT['BSL'], UP_PHOT['HMS'],\
        UP_PHOT['MCC_5'], UP_PHOT['TPR_5'], UP_PHOT['PPV_5'],\
        UP_PHOT['MCC_9'], UP_PHOT['TPR_9'], UP_PHOT['PPV_9']

    LSR_UPM, BSL_UPM, HMS_UPM, MCC_5_UPM, TPR_5_UPM, PPV_5_UPM, MCC_9_UPM,\
        TPR_9_UPM, PPV_9_UPM = UP_PM['LSR'], UP_PM['BSL'], UP_PM['HMS'],\
        UP_PM['MCC_5'], UP_PM['TPR_5'], UP_PM['PPV_5'], UP_PM['MCC_9'],\
        UP_PM['TPR_9'], UP_PM['PPV_9']

    pyU_PHOT = np.array([
        LSR_pUPH, BSL_pUPH, HMS_pUPH, MCC_5_pUPH, TPR_5_pUPH, PPV_5_pUPH,
        MCC_9_pUPH, TPR_9_pUPH, PPV_9_pUPH])
    pyU_PM = np.array([
        LSR_pUPM, BSL_pUPM, HMS_pUPM, MCC_5_pUPM, TPR_5_pUPM, PPV_5_pUPM,
        MCC_9_pUPM, TPR_9_pUPM, PPV_9_pUPM])
    U_PHOT = np.array([
        LSR_UPH, BSL_UPH, HMS_UPH, MCC_5_UPH, TPR_5_UPH, PPV_5_UPH, MCC_9_UPH,
        TPR_9_UPH, PPV_9_UPH])
    U_PM = np.array([
        LSR_UPM, BSL_UPM, HMS_UPM, MCC_5_UPM, TPR_5_UPM, PPV_5_UPM, MCC_9_UPM,
        TPR_9_UPM, PPV_9_UPM])

    comp_PHOT = pyU_PHOT - U_PHOT
    comp_PM = pyU_PM - U_PM

    win_PHOT, loss_PHOT, emp_PHOT, win_PM, loss_PM, emp_PM = np.zeros(9),\
        np.zeros(9), np.zeros(9), np.zeros(9), np.zeros(9), np.zeros(9)

    for i in range(len(pyU_PHOT)):
        for j in range(len(comp_PHOT[i])):
            if comp_PHOT[i][j] > tie_max:
                win_PHOT[i] = win_PHOT[i] + 1
            elif comp_PHOT[i][j] < tie_min:
                loss_PHOT[i] = loss_PHOT[i] + 1
            else:
                emp_PHOT[i] = emp_PHOT[i] + 1

    for i in range(len(pyU_PM)):
        for j in range(len(comp_PM[i])):
            if comp_PM[i][j] > tie_max:
                win_PM[i] = win_PM[i] + 1
            elif comp_PM[i][j] < tie_min:
                loss_PM[i] = loss_PM[i] + 1
            else:
                emp_PM[i] = emp_PM[i] + 1

    return CI_PM, CI_PHOT, comp_PHOT, comp_PM, win_PHOT, loss_PHOT, emp_PHOT,\
        win_PM, loss_PM, emp_PM


def makePlot(
    fold, tie_max, tie_min, H, m, CI_PM, CI_PHOT, comp_PHOT, comp_PM, win_PHOT,
        loss_PHOT, emp_PHOT, win_PM, loss_PM, emp_PM):
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
    file_out = fold + 'plots/H{}/'.format(H) + '{}.png'.format(m)
    plt.savefig(file_out, dpi=150, bbox_inches='tight')

    fig = plt.figure(figsize=(20, 10))
    CIPlot(tie_min, tie_max, CI_PM, comp_PM, CI_PHOT, comp_PHOT)
    file_out = fold + 'plots/CI_H{}/'.format(H) + '{}.png'.format(m)
    fig.tight_layout()
    plt.savefig(file_out, dpi=150, bbox_inches='tight')


def CIPlot(tie_min, tie_max, CI_PM, comp_PM, CI_PHOT, comp_PHOT):
    """
    """
    metrics = (
        'LSR', 'BSL', 'HMS', 'MCC_5', 'TPR_5', 'PPV_5', 'MCC_9', 'TPR_9',
        'PPV_9')

    for i, met in enumerate(metrics):
        ax = plt.subplot(int("33" + str(i + 1)))
        plt.title(met)
        plt.axhline(0., c='k', ls=':')
        plt.axhline(tie_min, ls=':', c='yellow', lw=1.5, zorder=1)
        plt.axhline(tie_max, ls=':', c='yellow', lw=1.5, zorder=1)
        plt.scatter(CI_PM, comp_PM[i], c='r', alpha=.5, label="PM", zorder=3)
        plt.scatter(
            CI_PHOT, comp_PHOT[i], c='b', alpha=.5, label="PHOT", zorder=3)
        if met in ('LSR', 'MCC_5', 'MCC_9'):
            plt.ylabel(r"$\Delta=(pyU - U)$")
        if met in ('MCC_9', 'TPR_9', 'PPV_9'):
            plt.xlabel("CI (log)")
        ax.set_xscale('log')
        if i == 0:
            plt.legend()


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
        text_color = 'white' if r * g * b < 0.5 else 'darkgrey'
        for y, (x, c) in enumerate(zip(xcenters, widths)):
            ax.text(x, y, str(int(c)), ha='center', va='center',
                    color=text_color)
    ax.legend(ncol=len(category_names), bbox_to_anchor=(0, 1),
              loc='lower left', fontsize='small')


def summaryPlot(fold, H, winloss_rates):
    """
    Summary of the combined metrics for all the methods.
    """
    fig = plt.figure(figsize=(15, 10))

    min_y, max_y = np.inf, 0
    plt.subplot(211)
    i = 0
    for k, v in winloss_rates.items():
        yoff = -.2 if (i % 2) == 0 else .17
        i += 1
        xr, yr = np.random.uniform(.015, .02, 2)
        plt.scatter(v[0] + xr, v[1] + yr, alpha=.7, s=50)
        # texts.append(plt.text(v[0], v[1], k, ha='center', va='center'))
        # plt.text(v[0], v[1], k, ha='center', va='center')
        plt.annotate(k, (v[0] - .025, v[1] + yoff))
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
    plt.ylim(0, 11)
    plt.legend()

    file_out = fold + 'plots/summary_{}.png'.format(H)
    fig.tight_layout()
    plt.savefig(file_out, dpi=150, bbox_inches='tight')


def readFiles(ruta=Path.cwd()):
    """
    Read files from the input folder
    """
    return [arch.name for arch in Path('input').iterdir() if arch.is_file()]


if __name__ == '__main__':
    main()
