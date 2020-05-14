
import configparser
from distutils.util import strtobool


def main():
    """
    """

    def vtype(var):
        tp, v = var.split('_')
        if tp == 'int':
            return int(v)
        elif tp == 'float':
            return float(v)
        elif tp == 'bool':
            return bool(strtobool(v))
        elif tp == 'str':
            return v

    # Read .ini file
    in_params = configparser.ConfigParser()
    in_params.read('params.ini')

    # Data columns
    ID_c, x_c, y_c = in_params['Data columns']['ID_coords'].split()
    data_cols = in_params['Data columns']['data'].split()
    data_errs = in_params['Data columns']['uncert'].split()

    # Arguments for the Outer Loop
    verbose, OL_runs, resampleFlag, PCAflag, PCAdims, prob_cnvrg =\
        vtype(in_params['Outer loop']['verbose']),\
        vtype(in_params['Outer loop']['OL_runs']),\
        vtype(in_params['Outer loop']['resampleFlag']),\
        vtype(in_params['Outer loop']['PCAflag']),\
        vtype(in_params['Outer loop']['PCAdims']),\
        vtype(in_params['Outer loop']['prob_cnvrg'])

    # Arguments for the Inner Loop
    clust_params = {
        'N_membs': vtype(in_params['Inner loop']['N_membs']),
        'N_cl_max': vtype(in_params['Inner loop']['N_cl_max'])}

    clust_method, C_thresh, RK_rad, RK_mode =\
        vtype(in_params['Inner loop']['clust_method']),\
        vtype(in_params['Inner loop']['C_thresh']),\
        vtype(in_params['Inner loop']['RK_rad']),\
        vtype(in_params['Inner loop']['RK_mode'])

    if RK_rad not in (.1, .3, .5, .7, .9):
        raise ValueError("RK radius must be one of the accepted "
                         "values: (.1, .2, .3, .4, .5, .6, .7, .8, .9)")

    cl_method_pars = {}
    for key, val in in_params['Clustering parameters'].items():
        cl_method_pars[key] = vtype(val)

    return ID_c, x_c, y_c, data_cols, data_errs, verbose, OL_runs,\
        resampleFlag, PCAflag, PCAdims, prob_cnvrg, clust_method,\
        C_thresh, RK_rad, RK_mode, clust_params, cl_method_pars
