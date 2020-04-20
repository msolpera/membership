
from modules.dataRead import read_data
from modules import member_index
# import modules.Voronoi_v4 as V5
import modules.Voronoi_v5 as V5
import modules.random_memb as RND
import modules.UP as UPMASK
from pathlib import Path


def main():
    """
    """
    methods = ('Voronoi_v5',) # 'random_memb', ('UPMASK',)

    # Create final table file
    with open("final_table.dat", "w") as f_out:
        f_out.write(
            'N_m,N_f,CI,CI_V,CI_BV,CI_pmRA,CI_pmDE,' +
            'C_V,P_V,L_V,C_R,P_R L_R,C_U,P_U,L_U\n')

    # Column names
    # For our synthetic clusters
    data_cols = ('pmRA', 'pmDE') #  , 'V', 'BV')
    data_errs = ('BV', 'BV') #  , 'BV', 'BV')
    # For UPMASK's synthetic clusters
    # data_cols = ('V', 'B_V', 'U_B', 'V_I')
    # data_errs = ('V', 'B_V', 'U_B', 'V_I')

    arch = readFiles('input')
    for file_name in arch:
        if file_name.startswith("README"):
            continue

        file_name = file_name.replace('up-', '')
        print("Processing: {}".format(file_name))

        CI = float(file_name[:4])
        CI_V = float(file_name[5:9])
        CI_BV = float(file_name[10:14])
        CI_pmRA = float(file_name[15:19])
        CI_pmDE = float(file_name[20:24])

        # Read data, normalized without outliers.
        ID, xy, data, data_err = read_data(
            'input/' + file_name, data_cols, data_errs)

        for met in methods:
            print('\nMethod:', met)

            if met == 'Voronoi_v5':
                memb_prob = V5.main(ID, xy, data, data_err)
                C, P, log_MI, Nm, Nf = member_index.main(ID, memb_prob)
                MI_v5 = [C, P, log_MI]

            if met == 'random_memb':
                memb_prob = RND.main(ID)
                C, P, log_MI = member_index.main(ID, memb_prob)[:3]
                print("random:", C, P, log_MI)
                MI_RD = [C, P, log_MI]

            if met == 'UPMASK':
                ID, memb_prob = UPMASK.main('input/UPMASK/up-' + file_name)
                C, P, log_MI = member_index.main(ID, memb_prob)[:3]
                print(C, P, log_MI)
                MI_UP = [C, P, log_MI]

        data_lst = (
            Nm, Nf, CI, CI_V, CI_BV, CI_pmRA, CI_pmDE, MI_v5[0], MI_v5[1],
            MI_v5[2], MI_RD[0], MI_RD[1], MI_RD[2], MI_UP[0], MI_UP[1],
            MI_UP[2])
        with open("final_table.dat", "a") as f_out:
            f_out.write(",".join(str(_) for _ in data_lst))
            f_out.write("\n")


def readFiles(ruta=Path.cwd()):
    """
    Read files from the input folder
    """
    return [arch.name for arch in Path('input').iterdir() if arch.is_file()]


if __name__ == '__main__':
    main()
