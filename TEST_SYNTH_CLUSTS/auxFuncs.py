
from astropy.table import Table
import numpy as np


# Folder where the files are located
fold = "metrics"

# Define the "tie" range.
tie_max = 0.005
tie_min = -1. * tie_max


def readTables(N_UPMASK, H, m, UP_alg=""):
    """
    It is important to use .group_by() to make sure that all the Tables are
    in the correct order
    """
    pyUP_PHOT = Table.read(
        fold + '/metrics_PHOT_{}_H_{}.dat'.format(m, H),
        format='ascii').group_by('Name')
    pyUP_PM = Table.read(
        fold + '/metrics_PM_{}_H_{}.dat'.format(m, H),
        format='ascii').group_by('Name')

    UP_PHOT = Table.read(
        fold + '/metrics_PHOT_UPMASK_600_' + UP_alg + str(N_UPMASK)
        + '_H_' + H + '.dat', format='ascii').group_by('Name')
    UP_PM = Table.read(
        fold + '/metrics_PM_UPMASK_600_' + UP_alg + str(N_UPMASK)
        + '_H_' + H + '.dat', format='ascii').group_by('Name')

    return pyUP_PHOT, pyUP_PM, UP_PHOT, UP_PM


def WinTieLoss(pyUP_PHOT, pyUP_PM, UP_PHOT, UP_PM, caller):
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

    if caller == 'metrics_bars':
        return win_PHOT, loss_PHOT, emp_PHOT, win_PM, loss_PM, emp_PM
    elif caller == 'summary':
        return CI_PM, CI_PHOT, win_PHOT, loss_PHOT, win_PM, loss_PM
    elif caller == 'CI':
        return CI_PM, CI_PHOT, comp_PHOT, comp_PM, pyU_PHOT, pyU_PM
