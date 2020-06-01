
import numpy as np
import pyUPMASK as V5
from pathlib import Path


def main():
    """
    """
    # # Print metrics for a given file in a selected folder
    # folder = 'output'
    # arch = readFiles(folder)
    # from astropy.io import ascii
    # data = ascii.read(folder + '/' + arch[1])
    # member_index(data['ID'], data['probability'])

    methods = ('pyUPMASK', 'UPMASK') # 'random_memb',  ('UPMASK',)

    # Create final table file
    with open("final_table.dat", "w") as f_out:
        f_out.write(
            'N_m,N_f,CI,CI_V,CI_BV,CI_pmRA,CI_pmDE,' +
            'C_V,P_V,L_V,C_R,P_R L_R,C_U,P_U,L_U\n')

    arch = readFiles('input')
    for file_name in arch:
        if file_name.startswith("README"):
            continue

        file_name = file_name.replace('up-', '')
        print("\n\nProcessing: {}".format(file_name))

        CI = float(file_name[:4])
        CI_V = float(file_name[5:9])
        CI_BV = float(file_name[10:14])
        CI_pmRA = float(file_name[15:19])
        CI_pmDE = float(file_name[20:24])

        for met in methods:

            if met == 'pyUPMASK':
                memb_prob = V5.main()
                print("\npyUPMASK")
                C, P, M, L, Nm, Nf = member_index(ID, memb_prob)
                MI_v5 = [C, P, M, L]

            if met == 'random_memb':
                memb_prob = RND.main(ID)
                C, P, M, L = member_index(ID, memb_prob)[:4]
                MI_RD = [C, P, M, L]

            if met == 'UPMASK':
                ID, memb_prob = UPMASK.main('input/UPMASK/up-' + file_name)
                print("\nUPMASK")
                C, P, M, L = member_index(ID, memb_prob)[:4]
                MI_UP = [C, P, M, L]

        # data_lst = (
        #     Nm, Nf, CI, CI_V, CI_BV, CI_pmRA, CI_pmDE, MI_v5[0], MI_v5[1],
        #     MI_v5[2], MI_RD[0], MI_RD[1], MI_RD[2], MI_UP[0], MI_UP[1],
        #     MI_UP[2])
        # with open("final_table.dat", "a") as f_out:
        #     f_out.write(",".join(str(_) for _ in data_lst))
        #     f_out.write("\n")


def readFiles(folder='input'):
    """
    Read files from the given folder
    """
    return [arch.name for arch in Path(folder).iterdir() if arch.is_file()]


def member_index(ID, memb_prob, MP_cut=(.5, .9)):
    """
    Obtain the Completeness (C), Purity (P), and Misclassification (M) using
    a cut on MP of 'MP_cut'.
    Also obtain the Logarithmic Scoring Rule (log_MI) which requires no cut.
    """
    # Avoid numerical errors with 0. and 1.
    memb_prob = np.clip(memb_prob, a_min=0.01, a_max=0.99)

    for MP in MP_cut:

        log_MI = 0.
        N_true, N_missed, N_selected, N_interlopers, N_field =\
            0., 0., 0., 0., 0.
        for i in range(len(memb_prob)):
            if str(ID[i])[0] == '1':
                N_true += 1.
                log_MI += np.log(memb_prob[i])
                if memb_prob[i] >= MP:
                    N_selected += 1
                else:
                    N_missed += 1
            else:
                N_field += 1
                log_MI += np.log(1. - memb_prob[i])
                if memb_prob[i] >= MP:
                    N_interlopers += 1

        # Make it so that 1 is the maximum (best) value, and shorten the range
        # with the logarithm
        L = 1. - np.log(1. - log_MI)

        # Completeness
        C = (N_true - N_missed) / N_true  # = N_selected / N_true
        # Purity
        if N_selected > 0:
            P = (N_selected - N_interlopers) / N_selected
        else:
            P = -np.inf
        # Misclassification
        M = 1. - (N_missed + N_interlopers) / N_true

        print("MP_cut={:.2f} -> C={:.3f}, P={:.3f}, M={:.3f}, L={:.3f}".format(
            MP, C, P, M, L))

    return C, P, M, L, N_true, N_field


def RND(ID):
    memb_prob = np.random.uniform(0., 1., len(ID))
    return memb_prob


def UPMASK(file_name):
    from astropy.table import Table
    data = Table.read(file_name, format='ascii')
    return data['ID'], data['probability']


if __name__ == '__main__':
    main()
