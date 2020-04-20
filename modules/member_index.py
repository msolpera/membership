
import numpy as np


def main(ID, memb_prob, MP_cut=.9):
    """
    Obtain the Completeness (C), Purity (P) using a cut on MP of 'MP_cut'.
    Also obtain the Logarithmic Scoring Rule (log_MI) which requires no cut.
    """
    # Avoid numerical errors with 0. and 1.
    memb_prob = np.clip(memb_prob, a_min=0.01, a_max=0.99)

    log_MI = 0.
    N_true, N_missed, N_selected, N_interlopers, N_field = 0., 0., 0., 0., 0.
    for i in range(len(memb_prob)):
        if str(ID[i])[0] == '1':
            N_true += 1.
            log_MI += np.log(memb_prob[i])
            if memb_prob[i] >= MP_cut:
                N_selected += 1
            else:
                N_missed += 1
        else:
            N_field += 1
            log_MI += np.log(1. - memb_prob[i])
            if memb_prob[i] >= MP_cut:
                N_interlopers += 1

    # Completeness
    C = (N_true - N_missed) / N_true
    # Purity
    if N_selected > 0:
        P = (N_selected - N_interlopers) / N_selected
    else:
        P = -np.inf

    return C, P, log_MI, N_true, N_field
