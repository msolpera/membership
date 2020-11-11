
import numpy as np
from pathlib import Path
from astropy.io import ascii
from astropy.table import Column
from astropy.table import Table
import configparser
from distutils.util import strtobool
from sklearn.preprocessing import MinMaxScaler
from .outlierRjct import stdRegion, sklearnMethod


def readINI():
    """
    Read .ini config file
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

    in_params = configparser.ConfigParser()
    in_params.read('params.ini')

    # Data columns
    data_columns = in_params['Data columns']
    ID_c = data_columns['ID']
    x_c, y_c = data_columns['xy_coords'].split()
    data_cols = data_columns['data'].split()
    oultr_method = data_columns.get('oultr_method')
    stdRegion_nstd = data_columns.getfloat('stdRegion_nstd')

    # Arguments for the Outer Loop
    outer_loop = in_params['Outer loop']
    rnd_seed, verbose, OL_runs, parallel_flag, parallel_procs, resampleFlag,\
        PCAflag, PCAdims, GUMM_flag, KDEP_flag =\
        outer_loop.get('rnd_seed'), outer_loop.getint('verbose'),\
        outer_loop.getint('OL_runs'), outer_loop.getboolean('parallel'),\
        outer_loop.get('processes'), outer_loop.getboolean('resampleFlag'),\
        outer_loop.getboolean('PCAflag'), outer_loop.getint('PCAdims'),\
        outer_loop.getboolean('GUMM_flag'), outer_loop.getboolean('KDEP_flag')
    GUMM_perc = outer_loop.get('GUMM_perc')
    if GUMM_perc != 'auto':
        GUMM_perc = float(GUMM_perc)

    # Only read if the code is set to re-sample the data.
    data_errs = []
    if resampleFlag:
        data_errs = data_columns['uncert'].split()

    # Arguments for the Inner Loop
    inner_loop = in_params['Inner loop']
    N_membs, clust_method, clRjctMethod, C_thresh =\
        inner_loop.getint('N_membs'), inner_loop.get('clust_method'),\
        inner_loop.get('clRjctMethod'), inner_loop.getfloat('C_thresh')

    if clRjctMethod not in ('rkfunc', 'kdetest', 'kdetestpy'):
        raise ValueError("'{}' is not a valid choice for clRjctMethod".format(
            clRjctMethod))

    cl_method_pars = {}
    for key, val in in_params['Clustering parameters'].items():
        cl_method_pars[key] = vtype(val)

    return ID_c, x_c, y_c, data_cols, data_errs, oultr_method, stdRegion_nstd,\
        rnd_seed, verbose, OL_runs, parallel_flag, parallel_procs,\
        resampleFlag, PCAflag, PCAdims, GUMM_flag, GUMM_perc, KDEP_flag,\
        N_membs, clust_method, clRjctMethod, C_thresh, cl_method_pars


def dread(file_path, ID_c, x_c, y_c, data_cols, data_errs):
    """
    """

    data = Table.read(file_path, format='ascii')
    N_d = len(data)
    print("Stars read         : {}".format(N_d))

    # Remove stars with no valid data
    try:
        msk = np.logical_or.reduce([~data[_].mask for _ in data_cols])
        data_rjct = data[~msk]
        data = data[msk]
        print("Stars removed      : {}".format(N_d - len(data)))
    except AttributeError:
        # No masked columns
        data_rjct = []
        pass

    # Separate data into groups
    if ID_c == 'None':
        N_d = len(data)
        ID_data = np.arange(1, N_d + 1)
    else:
        ID_data = data[ID_c]
    xy_data, cl_data = np.array([data[x_c], data[y_c]]).T,\
        np.array([data[_] for _ in data_cols]).T

    cl_errs = np.array([])
    if data_errs:
        cl_errs = np.array([data[_] for _ in data_errs]).T
    print("Data dimensions    : {}".format(cl_data.shape[1]))

    return data, ID_data, xy_data, cl_data, cl_errs, data_rjct


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


def dwrite(out_folder, file_path, full_data, msk_data, probs_all, probs_mean):
    """
    """
    if msk_data is not None:
        out_path = Path(out_folder, *file_path.parts[1:])

        for i, p in enumerate(probs_all):
            # Fill masked data with '-1'
            p0 = np.zeros(len(full_data)) - 1.
            p0[msk_data] = p
            full_data.add_column(Column(np.round(p0, 2), name='prob' + str(i)))

        pf = np.zeros(len(full_data)) - 1.
        pf[msk_data] = probs_mean
        full_data.add_column(Column(np.round(pf, 2), name='probs_final'))
    else:
        ext = file_path.suffix
        file_path = Path(str(file_path).replace(ext, '_rjct' + ext))
        out_path = Path(out_folder, *file_path.parts[1:])

    # Create sub-folder if it does not exist
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ascii.write(full_data, out_path, overwrite=True)


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
