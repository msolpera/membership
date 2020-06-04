
from astropy.io import ascii
import numpy as np
from sklearn.metrics import log_loss, brier_score_loss, roc_auc_score,\
    matthews_corrcoef, recall_score, precision_score


def main():
    """
    """
    # List of clusters already process with pyUPMASK or UPMASK
    clusts_lst = ["./aux_funcs/UPMASK/up-RESULTS.dat"]

    final_dct = {
        'Name': [], 'CI': [], 'LSR': [], 'BSL': [], 'AUC': [],
        'MCC_5': [], 'TPR_5': [], 'PPV_5': [],
        'MCC_9': [], 'TPR_9': [], 'PPV_9': []}

    for fpath in clusts_lst:
        data = ascii.read(fpath)
        metrics_dct = getMetrics(data, verbose=False)

        fname = fpath.split('/')[-1].split('.')[0]
        final_dct['Name'].append(fname)
        for k, v in metrics_dct.items():
            final_dct[k].append(metrics_dct[k])

    ascii.write(final_dct, "metrics.dat", format='csv', overwrite=True)


def getMetrics(data, eps=1e-2, verbose=True):
    """
    Obtain performance metrics that require no probability cut, as well
    as others that do require one.

    UPMASK's metrics are badly defined. their TPR is in fact PPV, and the
    MMR is the (proper) TPR
    """

    try:
        # For our synthetic clusters
        ID = data['ID']
        memb_id = '1'
    except KeyError:
        # For UPMASK synthetic clusters
        ID = data['field']
        memb_id = '0'

    try:
        # Processed with pyUPMASK
        memb_prob = data['probs_final']
        # pyUPMASK assigns P=-1 to outliers
        memb_prob = np.clip(memb_prob, a_min=0., a_max=1.)
    except KeyError:
        # Processed with UPMASK
        memb_prob = data['probability']

    # Members are identified with a '1'
    y_true = []
    for ID_v in ID:
        if str(ID_v).startswith(memb_id):
            y_true.append(1)
        else:
            y_true.append(0)

    N_membs = sum(y_true)
    N_field = len(ID) - N_membs
    if verbose:
        print("N_field/N_membs={:.3f}".format(N_field / N_membs))

    # Invert so that 1 is the maximum (best) value. Avoid numerical errors
    # with 0. and 1 by clipping
    LSR = 1. - log_loss(y_true, memb_prob, eps=eps)
    # Brier score loss. Invert so that 1 is the maximum (best) value.
    BSL = 1. - brier_score_loss(y_true, memb_prob)
    # Area Under the Receiver Operating Characteristic Curve (ROC AUC)
    AUC = roc_auc_score(y_true, memb_prob)

    if verbose:
        print("LSR={:.3f}, BSL={:.3f}, AUC={:.3f}".format(LSR, BSL, AUC))

    # Dictionary with metrics
    metrics_dct = {'CI': N_field / N_membs, 'LSR': LSR, 'BSL': BSL, 'AUC': AUC}

    # Calculate the metrics below for the following probability cuts.
    Prob_cuts = (.5, .9)

    # True Positive  (TP) : correctly identified
    # False Positive (FP) : incorrectly identified; Type 1 Error
    # True Negative  (TN) : correctly rejected
    # False Negative (FN) : incorrectly rejected; Type 2 Error
    #
    # N_true = TP + FN
    # N_field = FP + TN

    for MP_cut in Prob_cuts:
        # FN, TP, FP, TN = 0., 0., 0., 0.
        y_pred = []
        for i, MP in enumerate(memb_prob):

            # This is a real member
            if y_true[i] == 1:
                if MP >= MP_cut:
                    # TP += 1
                    y_pred.append(1)
                else:
                    # FN += 1
                    y_pred.append(0)

            # This is a field star
            else:
                if MP >= MP_cut:
                    # FP += 1
                    y_pred.append(1)
                else:
                    # TN += 1
                    y_pred.append(0)

        # Matthews correlation coefficient
        MCC = matthews_corrcoef(y_true, y_pred)
        # Sensitivity (recall, hit rate, true positive rate)
        TPR = recall_score(y_true, y_pred)
        # Precision (positive predictive value)
        PPV = precision_score(y_true, y_pred)

        metrics_dct['MCC_{}'.format(int(MP_cut * 10))] = MCC
        metrics_dct['TPR_{}'.format(int(MP_cut * 10))] = TPR
        metrics_dct['PPV_{}'.format(int(MP_cut * 10))] = PPV

        if verbose:
            print("P_c={}, MCC={:.3f}, TPR={:.3f}, PPV={:.3f}".format(
                MP_cut, MCC, TPR, PPV))

    return metrics_dct


if __name__ == '__main__':
    main()
