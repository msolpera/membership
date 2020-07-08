
from pathlib import Path
from astropy.io import ascii
import numpy as np
from sklearn.metrics import log_loss, brier_score_loss,\
    matthews_corrcoef, recall_score, precision_score


def main(verbose=False):
    """
    Obtain performance metrics for the clusters processed with pyUPMASK or
    UPMASK.

    The input files (which are the outputs from pyUPMASk or UPMASK) are read
    from the 'output/' folder.
    """
    # mode = ('UPMASK_600', 'pyUPMASK_600')
    mode = ('autoperc_inner_GUMM7',)
    features = ('PHOT', 'PM')
    Hval = ('auto',)  # 'symm', 'SR05')

    for H in Hval:
        print("\n{}".format(H))
        for m in mode:
            print(m)
            for f in features:
                print(f)
                subfold = m + '/' + f

                final_dct = {
                    'Name': [], 'CI': [], 'LSR': [], 'BSL': [], 'HMS': [],
                    'MCC_5': [], 'TPR_5': [], 'PPV_5': [],
                    'MCC_9': [], 'TPR_9': [], 'PPV_9': []}

                outp = Path('../output/' + subfold)
                for fpath in outp.iterdir():
                    fname = '_'.join(fpath.name.split('/')[-1].split('.')[:-1])
                    if fname == '':
                        continue
                    print("  ", fname)

                    data = ascii.read(fpath)
                    metrics_dct = getMetrics(data, H, verbose)

                    final_dct['Name'].append(fname)
                    for k, v in metrics_dct.items():
                        final_dct[k].append(metrics_dct[k])

                outf = Path('../TEST_SYNTH_CLUSTS/test_results/').joinpath(
                    'metrics_{}_{}_H_{}.dat'.format(f, m, H))
                ascii.write(final_dct, outf, format='csv', overwrite=True)


def getMetrics(data, Hval, verbose, eps=1e-2):
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
        # memb_prob = data['prob0']
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

    # Invert so that 1 is the maximum (best) value. The 'eps' parameter is
    # there to  prevent numerical errors with the logarithm of 0. or 1.
    LSR = 1. - log_loss(y_true, memb_prob, eps=eps)
    # Brier score loss. Invert so that 1 is the maximum (best) value.
    BSL = 1. - brier_score_loss(y_true, memb_prob)
    # Area Under the Receiver Operating Characteristic Curve (ROC AUC)
    # AUC = roc_auc_score(y_true, memb_prob)

    from rpy2.robjects import r
    from rpy2.robjects import numpy2ri
    from rpy2.robjects.packages import importr
    importr('hmeasure')
    numpy2ri.activate()
    r.assign('true_labels', np.array(y_true))
    r.assign('scores', np.array([memb_prob.data, ]).T)
    if Hval == 'auto':
        # Auto severity.ratio (=N_1/N_0)
        r('results <- HMeasure(true_labels, scores)')
    elif Hval == 'symm':
        # Symmetric severity.ratio
        r('results <- HMeasure(true_labels, scores, severity.ratio=1)')
    elif Hval == 'SR05':
        r('results <- HMeasure(true_labels, scores, severity.ratio=0.5)')
    HMS = list(r('results$metrics'))[0][0]

    if verbose:
        print("LSR={:.3f}, BSL={:.3f}, HMS={:.3f}".format(LSR, BSL, HMS))

    # Dictionary with metrics
    metrics_dct = {'CI': N_field / N_membs, 'LSR': LSR, 'BSL': BSL, 'HMS': HMS}

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
