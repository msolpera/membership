import matplotlib.pyplot as plt
from auxFuncs import tie_min, tie_max, readTables
import numpy as np

fold = '/metrics'


def main():
    """
    Plot the methods versus CI for all the clusters.
    """
    configs = ("Voron", "kNNde", "Agglo", 'MiniB', "KMean", "Gauss")
    Hval = 'auto'  # 'symm', 'SR05'
    N_UPMASK = "25"  # "50"

    CI_PM, CI_PHOT, comp_PHOT, comp_PM, pyU_PHOT, pyU_PM =\
        [], [], [], [], [], []
    for i, m in enumerate(configs):
        print(" ", m)
        pyUP_PHOT, pyUP_PM, UP_PHOT, UP_PM = readTables(N_UPMASK, Hval, m)

        values = WinTieLoss(
            tie_max, tie_min, pyUP_PHOT, pyUP_PM, UP_PHOT, UP_PM, 'CI')
        CI_PM.append(values[0])
        CI_PHOT.append(values[1])
        comp_PHOT.append(values[2])
        comp_PM.append(values[3])
        pyU_PHOT.append(values[4])
        pyU_PM.append(values[5])
    conf = ("VOR", "KNN", "AGG", 'MBK', "KMS", "GMM")
    fig = plt.figure(figsize=(13, 10))
    sct_rdm_pm = np.random.normal(-2, 2, len(CI_PM[0]))
    sct_rdm_ph = np.random.normal(-2, 2, len(CI_PHOT[0]))
    for j, met in enumerate(conf):
        ax = plt.subplot(3, 2, j + 1)
        plt.title(met)
        plt.axhline(0., c='k', ls=':')
        x_pm = abs(CI_PM[j] + sct_rdm_pm)
        x_ph = abs(CI_PHOT[j] + sct_rdm_ph)
        plt.scatter(
            x_pm, comp_PM[j], c='darkgray', alpha=.3, label="PM", zorder=3,
            edgecolors='black')
        plt.scatter(
            x_ph, comp_PHOT[j], c='b', marker='^', alpha=.3, label="PHOT",
            zorder=3, edgecolors='darkblue')
        ymin, ymax = ax.get_ylim()

        ax.axhspan(tie_min, tie_max, alpha=.1, color='yellow')
        ax.axhspan(-10, tie_min, alpha=.1, color='red')
        ax.axhspan(tie_max, 10, alpha=.1, color='green')
        plt.ylim(ymin, ymax)

        if met in ('VOR', 'AGG', 'KMS'):
            plt.ylabel("pyUPMASK - UPMASK")
        if met in ('KMS', 'GMM'):
            plt.xlabel("log(CI)")
        ax.set_xscale('log')
        if j == 0:
            plt.legend()
        fig.tight_layout()
    # plt.savefig('CI_pyU.png', dpi=300, bbox_inches='tight')
    plt.show()


def WinTieLoss(tie_max, tie_min, pyUP_PHOT, pyUP_PM, UP_PHOT, UP_PM, caller):
    """
    """
    CI_PM, CI_PHOT = np.repeat(pyUP_PM['CI'], 9), np.repeat(pyUP_PHOT['CI'], 9)

    LSR_pUPH, BSL_pUPH, HMS_pUPH, MCC_5_pUPH, TPR_5_pUPH, PPV_5_pUPH,\
        MCC_9_pUPH, TPR_9_pUPH, PPV_9_pUPH = list(pyUP_PHOT['LSR']),\
        list(pyUP_PHOT['BSL']), list(pyUP_PHOT['HMS']), list(pyUP_PHOT['MCC_5']),\
        list(pyUP_PHOT['TPR_5']), list(pyUP_PHOT['PPV_5']), list(pyUP_PHOT['MCC_9']),\
        list(pyUP_PHOT['TPR_9']), list(pyUP_PHOT['PPV_9'])

    LSR_pUPM, BSL_pUPM, HMS_pUPM, MCC_5_pUPM, TPR_5_pUPM, PPV_5_pUPM,\
        MCC_9_pUPM, TPR_9_pUPM, PPV_9_pUPM = list(pyUP_PM['LSR']),\
        list(pyUP_PM['BSL']), list(pyUP_PM['HMS']), list(pyUP_PM['MCC_5']),\
        list(pyUP_PM['TPR_5']), list(pyUP_PM['PPV_5']), list(pyUP_PM['MCC_9']),\
        list(pyUP_PM['TPR_9']), list(pyUP_PM['PPV_9'])

    LSR_UPH, BSL_UPH, HMS_UPH, MCC_5_UPH, TPR_5_UPH, PPV_5_UPH, MCC_9_UPH,\
        TPR_9_UPH, PPV_9_UPH = list(UP_PHOT['LSR']),\
        list(UP_PHOT['BSL']), list(UP_PHOT['HMS']), list(UP_PHOT['MCC_5']),\
        list(UP_PHOT['TPR_5']), list(UP_PHOT['PPV_5']), list(UP_PHOT['MCC_9']),\
        list(UP_PHOT['TPR_9']), list(UP_PHOT['PPV_9'])

    LSR_UPM, BSL_UPM, HMS_UPM, MCC_5_UPM, TPR_5_UPM, PPV_5_UPM, MCC_9_UPM,\
        TPR_9_UPM, PPV_9_UPM = list(UP_PM['LSR']),\
        list(UP_PM['BSL']), list(UP_PM['HMS']), list(UP_PM['MCC_5']),\
        list(UP_PM['TPR_5']), list(UP_PM['PPV_5']), list(UP_PM['MCC_9']),\
        list(UP_PM['TPR_9']), list(UP_PM['PPV_9'])

    pyU_PHOT = LSR_pUPH + BSL_pUPH + HMS_pUPH + MCC_5_pUPH + TPR_5_pUPH +\
        PPV_5_pUPH + MCC_9_pUPH + TPR_9_pUPH + PPV_9_pUPH
    pyU_PM = LSR_pUPM + BSL_pUPM + HMS_pUPM + MCC_5_pUPM + TPR_5_pUPM +\
        PPV_5_pUPM + MCC_9_pUPM + TPR_9_pUPM + PPV_9_pUPM
    U_PHOT = LSR_UPH + BSL_UPH + HMS_UPH + MCC_5_UPH + TPR_5_UPH +\
        PPV_5_UPH + MCC_9_UPH + TPR_9_UPH + PPV_9_UPH
    U_PM = LSR_UPM + BSL_UPM + HMS_UPM + MCC_5_UPM + TPR_5_UPM +\
        PPV_5_UPM + MCC_9_UPM + TPR_9_UPM + PPV_9_UPM

    comp_PHOT = np.array(pyU_PHOT) - np.array(U_PHOT)
    comp_PM = np.array(pyU_PM) - np.array(U_PM)

    win_PHOT, loss_PHOT, emp_PHOT, win_PM, loss_PM, emp_PM = np.zeros(1),\
        np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1)

    for i in range(len(comp_PHOT)):
        if comp_PHOT[i] > tie_max:
            win_PHOT = win_PHOT + 1
        elif comp_PHOT[i] < tie_min:
            loss_PHOT = loss_PHOT + 1
        else:
            emp_PHOT = emp_PHOT + 1

    for j in range(len(comp_PM)):
        if comp_PM[j] > tie_max:
            win_PM = win_PM + 1
        elif comp_PM[j] < tie_min:
            loss_PM = loss_PM + 1
        else:
            emp_PM = emp_PM + 1

    if caller == 'metrics_bars':
        return win_PHOT, loss_PHOT, emp_PHOT, win_PM, loss_PM, emp_PM
    elif caller == 'summary':
        return CI_PM, CI_PHOT, win_PHOT, loss_PHOT, win_PM, loss_PM
    elif caller == 'CI':
        return CI_PM, CI_PHOT, comp_PHOT, comp_PM, pyU_PHOT, pyU_PM


if __name__ == '__main__':
    main()
