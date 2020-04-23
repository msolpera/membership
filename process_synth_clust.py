
import configparser
from modules.dataRead import read_data
from modules import member_index
# import modules.Voronoi_v4 as V5
import modules.pyUPMASK as V5
import modules.random_memb as RND
import modules.UP as UPMASK
from pathlib import Path


def main():
    """
    """
    methods = ('pyUPMASK',) # 'random_memb', ('UPMASK',)

    # Create final table file
    with open("final_table.dat", "w") as f_out:
        f_out.write(
            'N_m,N_f,CI,CI_V,CI_BV,CI_pmRA,CI_pmDE,' +
            'C_V,P_V,L_V,C_R,P_R L_R,C_U,P_U,L_U\n')

    ID_c, x_c, y_c, data_cols, data_errs, OL_runs, resampleFlag, PCAflag,\
        PCAdims, otlrFlag, prob_cnvrg, clust_method, otlrFlag, C_thresh,\
        unif_method, RK_rad, clust_params = readINI()

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
            'input/' + file_name, ID_c, x_c, y_c, data_cols, data_errs)

        for met in methods:
            print('\nMethod:', met)

            if met == 'pyUPMASK':
                memb_prob = V5.main(
                    ID, xy, data, data_err, OL_runs, resampleFlag, PCAflag,
                    PCAdims, prob_cnvrg, clust_method, otlrFlag, C_thresh,
                    unif_method, RK_rad, clust_params)
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


def readINI():
    """
    """
    # Read .ini file
    in_params = configparser.ConfigParser()
    in_params.read('params.ini')

    # Data columns
    ID_c, x_c, y_c = in_params['Data columns']['ID_coords'].split()
    data_cols = in_params['Data columns']['data'].split()
    data_errs = in_params['Data columns']['uncert'].split()

    # Arguments for the Outer Loop
    OL_runs, resampleFlag, PCAflag, PCAdims, otlrFlag, prob_cnvrg =\
        int(in_params['Outer loop']['OL_runs']),\
        bool(in_params['Outer loop']['resampleFlag']),\
        bool(in_params['Outer loop']['PCAflag']),\
        in_params['Outer loop']['PCAdims'],\
        bool(in_params['Outer loop']['otlrFlag']),\
        float(in_params['Outer loop']['prob_cnvrg'])

    # Arguments for the Inner Loop
    clust_method, C_thresh, unif_method, RK_rad =\
        in_params['Inner loop']['clust_method'],\
        float(in_params['Inner loop']['C_thresh']),\
        in_params['Inner loop']['unif_method'],\
        float(in_params['Inner loop']['RK_rad'])

    clust_params = in_params['General clustering parameters']

    return ID_c, x_c, y_c, data_cols, data_errs,\
        OL_runs, resampleFlag, PCAflag, PCAdims, otlrFlag, prob_cnvrg,\
        clust_method, otlrFlag, C_thresh, unif_method, RK_rad, clust_params


def readFiles(ruta=Path.cwd()):
    """
    Read files from the input folder
    """
    return [arch.name for arch in Path('input').iterdir() if arch.is_file()]


if __name__ == '__main__':
    main()
