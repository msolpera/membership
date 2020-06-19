from pathlib import Path
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
# from matplotlib.pyplot import figure
# from astropy.io import ascii


def main():
    """
    Plot the statistics obtained after applying the getMetrics() script to
    the pyUPMASK and UPMASK results.
    """

    # Define the "tie" range.
    tie_max = .05
    tie_min = -1. * tie_max

    # Folder where the files are located
    fold = "../TEST_SYNTH_CLUSTS/test_results/"
    pyUP_PHOT = Table.read(fold + 'metrics_pyUP-PHOT.dat', format='ascii')
    pyUP_PM = Table.read(fold + 'metrics_pyUP-PM.dat', format='ascii')
    UP_PHOT = Table.read(fold + 'metrics_UP-PHOT.dat', format='ascii')
    UP_PM = Table.read(fold + 'metrics_UP-PM.dat', format='ascii')

    LSR_pUPH, BSL_pUPH, AUC_pUPH, MCC_5_pUPH, TPR_5_pUPH, PPV_5_pUPH,\
        MCC_9_pUPH, TPR_9_pUPH, PPV_9_pUPH = pyUP_PHOT['LSR'],\
        pyUP_PHOT['BSL'], pyUP_PHOT['AUC'], pyUP_PHOT['MCC_5'],\
        pyUP_PHOT['TPR_5'], pyUP_PHOT['PPV_5'], pyUP_PHOT['MCC_9'],\
        pyUP_PHOT['TPR_9'], pyUP_PHOT['PPV_9']

    LSR_pUPM, BSL_pUPM, AUC_pUPM, MCC_5_pUPM, TPR_5_pUPM, PPV_5_pUPM,\
        MCC_9_pUPM, TPR_9_pUPM, PPV_9_pUPM = pyUP_PM['LSR'], pyUP_PM['BSL'],\
        pyUP_PM['AUC'], pyUP_PM['MCC_5'], pyUP_PM['TPR_5'], pyUP_PM['PPV_5'],\
        pyUP_PM['MCC_9'], pyUP_PM['TPR_9'], pyUP_PM['PPV_9']

    LSR_UPH, BSL_UPH, AUC_UPH, MCC_5_UPH, TPR_5_UPH, PPV_5_UPH, MCC_9_UPH,\
        TPR_9_UPH, PPV_9_UPH = UP_PHOT['LSR'], UP_PHOT['BSL'], UP_PHOT['AUC'],\
        UP_PHOT['MCC_5'], UP_PHOT['TPR_5'], UP_PHOT['PPV_5'],\
        UP_PHOT['MCC_9'], UP_PHOT['TPR_9'], UP_PHOT['PPV_9']

    LSR_UPM, BSL_UPM, AUC_UPM, MCC_5_UPM, TPR_5_UPM, PPV_5_UPM, MCC_9_UPM,\
        TPR_9_UPM, PPV_9_UPM = UP_PM['LSR'], UP_PM['BSL'], UP_PM['AUC'],\
        UP_PM['MCC_5'], UP_PM['TPR_5'], UP_PM['PPV_5'], UP_PM['MCC_9'],\
        UP_PM['TPR_9'], UP_PM['PPV_9']

    pyU_PHOT = np.array([
        LSR_pUPH, BSL_pUPH, AUC_pUPH, MCC_5_pUPH, TPR_5_pUPH, PPV_5_pUPH,
        MCC_9_pUPH, TPR_9_pUPH, PPV_9_pUPH])
    pyU_PM = np.array([
        LSR_pUPM, BSL_pUPM, AUC_pUPM, MCC_5_pUPM, TPR_5_pUPM, PPV_5_pUPM,
        MCC_9_pUPM, TPR_9_pUPM, PPV_9_pUPM])
    U_PHOT = np.array([
        LSR_UPH, BSL_UPH, AUC_UPH, MCC_5_UPH, TPR_5_UPH, PPV_5_UPH, MCC_9_UPH,
        TPR_9_UPH, PPV_9_UPH])
    U_PM = np.array([
        LSR_UPM, BSL_UPM, AUC_UPM, MCC_5_UPM, TPR_5_UPM, PPV_5_UPM, MCC_9_UPM,
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

    # PLOTS
    category_names_PM = ['Loss_PM', 'Tie_PM', 'Win_PM']
    Results_PM = {
        'LSR': [loss_PM[0], emp_PM[0], win_PM[0]],
        'BSL': [loss_PM[1], emp_PM[1], win_PM[1]],
        'AUC': [loss_PM[2], emp_PM[2], win_PM[2]],
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
        'AUC': [loss_PHOT[2], emp_PHOT[2], win_PHOT[2]],
        'MCC_5': [loss_PHOT[3], emp_PHOT[3], win_PHOT[3]],
        'TPR_5': [loss_PHOT[4], emp_PHOT[4], win_PHOT[4]],
        'PPV_5': [loss_PHOT[5], emp_PHOT[5], win_PHOT[5]],
        'MCC_9': [loss_PHOT[6], emp_PHOT[6], win_PHOT[6]],
        'TPR_9': [loss_PHOT[7], emp_PHOT[7], win_PHOT[7]],
        'PPV_9': [loss_PHOT[8], emp_PHOT[8], win_PHOT[8]]
    }

    ax = plt.subplot(211)
    survey(ax, Results_PM, category_names_PM)
    ax = plt.subplot(212)
    survey(ax, Results_PHOT, category_names_PHOT)
    plt.show()


def survey(ax, results, category_names):
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

    # fig, ax = plt.subplots(figsize=(9.2, 3))
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

    # return fig, ax


def readFiles(ruta=Path.cwd()):
    """
    Read files from the input folder
    """
    return [arch.name for arch in Path('input').iterdir() if arch.is_file()]


if __name__ == '__main__':
    main()
