
import numpy as np
import pyUPMASK as V5
from pathlib import Path


def main():
    """
    """

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


def RND(ID):
    memb_prob = np.random.uniform(0., 1., len(ID))
    return memb_prob


def UPMASK(file_name):
    from astropy.table import Table
    data = Table.read(file_name, format='ascii')
    return data['ID'], data['probability']


if __name__ == '__main__':
    main()