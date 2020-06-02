
import numpy as np
from astropy.io import ascii
from astropy.table import Column
from astropy.table import Table
import configparser
from distutils.util import strtobool
from sklearn.preprocessing import MinMaxScaler
from .outlierRjct import stdRegion, sklearnMethod


def readINI():
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
    oultr_method = vtype(in_params['Data columns']['oultr_method'])
    stdRegion_nstd = vtype(in_params['Data columns']['stdRegion_nstd'])

    # Arguments for the Outer Loop
    rnd_seed, verbose, OL_runs, resampleFlag, PCAflag, PCAdims, GUMM_flag,\
        GUMM_perc =\
        vtype(in_params['Outer loop']['rnd_seed']),\
        vtype(in_params['Outer loop']['verbose']),\
        vtype(in_params['Outer loop']['OL_runs']),\
        vtype(in_params['Outer loop']['resampleFlag']),\
        vtype(in_params['Outer loop']['PCAflag']),\
        vtype(in_params['Outer loop']['PCAdims']),\
        vtype(in_params['Outer loop']['GUMM_flag']),\
        vtype(in_params['Outer loop']['GUMM_perc']),\

    # Only read if the code is set to re-sample the data.
    data_errs = []
    if resampleFlag:
        data_errs = in_params['Data columns']['uncert'].split()

    # Arguments for the Inner Loop
    N_membs, clust_method, clRjctMethod, RK_rad, C_thresh =\
        vtype(in_params['Inner loop']['N_membs']),\
        vtype(in_params['Inner loop']['clust_method']),\
        vtype(in_params['Inner loop']['clRjctMethod']),\
        vtype(in_params['Inner loop']['RK_rad']),\
        vtype(in_params['Inner loop']['C_thresh'])

    if clRjctMethod not in ('rkfunc', 'kdetest', 'kdetestpy'):
        raise ValueError("'{}' is not a valid choice for clRjctMethod".format(
            clRjctMethod))
    if RK_rad not in (.3, .5, .7):
        raise ValueError("RK radius must be one of the accepted "
                         "values: (.1, .2, .3, .4, .5, .6, .7, .8, .9)")

    cl_method_pars = {}
    for key, val in in_params['Clustering parameters'].items():
        cl_method_pars[key] = vtype(val)

    return ID_c, x_c, y_c, data_cols, data_errs, oultr_method, stdRegion_nstd,\
        rnd_seed, verbose, OL_runs, resampleFlag, PCAflag, PCAdims, GUMM_flag,\
        GUMM_perc, N_membs, clust_method, clRjctMethod, RK_rad, C_thresh,\
        cl_method_pars


def dread(
    file_name, ID_c, x_c, y_c, data_cols, data_errs, oultr_method,
        stdRegion_nstd):
    """
    """

    data = Table.read(file_name, format='ascii')
    print("\nStars read         : {}".format(len(data)))

    # Separate data into groups
    ID_data, xy_data, cl_data = data[ID_c],\
        np.array([data[x_c], data[y_c]]).T,\
        np.array([data[_] for _ in data_cols]).T

    cl_errs = np.array([])
    if data_errs:
        cl_errs = np.array([data[_] for _ in data_errs]).T

    return data, ID_data, xy_data, cl_data, cl_errs


def dmask(ID, xy, pdata, perrs, oultr_method, stdRegion_nstd):
    """
    """
    if oultr_method == 'stdregion':
        msk_data = stdRegion(pdata, stdRegion_nstd)
    else:
        msk_data = sklearnMethod(pdata, oultr_method)

    ID_data, xy_data, cl_data = ID[msk_data], xy[msk_data], pdata[msk_data]

    if perrs.any():
        data_err = perrs[msk_data]
    else:
        data_err = np.array([])

    print("Masked outliers    : {}".format((~msk_data).sum()))
    if oultr_method == 'stdregion':
        print(" N_std             : {}".format(stdRegion_nstd))

    return msk_data, ID_data, xy_data, cl_data, data_err


def dxynorm(xy_data):
    """
    """
    _xrange, _yrange = np.ptp(xy_data, 0)
    perc_sq = abs(1. - _xrange / _yrange)
    if perc_sq > .05:
        print((
            "WARNING: (x, y) frame deviates from a square region "
            "by {:.0f}%").format(perc_sq * 100.))
    xy = MinMaxScaler().fit(xy_data).transform(xy_data)
    print("Coordinates scaled : [0, 1]")

    return xy


def dwrite(file_name, full_data, msk_data, probs_all, probs_mean):
    """
    """
    fout = './output/' + file_name
    for i, p in enumerate(probs_all):
        # Fill masked data with '-1'
        p0 = np.zeros(len(full_data)) - 1.
        p0[msk_data] = p
        full_data.add_column(Column(np.round(p0, 2), name='prob' + str(i)))

    pf = np.zeros(len(full_data)) - 1.
    pf[msk_data] = probs_mean
    full_data.add_column(Column(np.round(pf, 2), name='probs_final'))

    ascii.write(full_data, fout, overwrite=True)


# def dataNorm(data_arr, err_data=None):
#     """
#     """
#     data_norm, err_norm = [], []
#     for i, arr in enumerate(data_arr.T):
#         min_array, max_array = np.nanmin(arr), np.nanmax(arr)
#         arr_delta = max_array - min_array
#         data_norm.append((arr - min_array) / arr_delta)

#         if err_data is not None:
#             err_norm.append(err_data.T[i] / arr_delta)

#         # # This normalization tends to make things more difficult
#         # mean_arr, std_arr = np.mean(arr[msk_data]), np.std(arr[msk_data])
#         # data_norm.append((arr[msk_data] - mean_arr) / std_arr)

#     return np.array(data_norm).T, np.array(err_norm).T
