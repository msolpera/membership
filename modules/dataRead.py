
import numpy as np
from astropy.table import Table
from sklearn.preprocessing import MinMaxScaler


def read_data(file_name, data_cols, data_errs):
    """
    """

    data = Table.read(file_name, format='ascii')
    print("\nN={} stars read".format(len(data)))
    # data = data['ID', 'x', 'y', 'V', 'BV', 'pmRA', 'pmDE']
    # data.remove_rows(np.where([c.data for c in data.mask.itercols()])[-1])
    # msk = data['Gmag'] < 16
    # data = data[msk]

    # Separate data into groups
    ID_data, xy_data, cl_data, data_err = data['ID'],\
        np.array([data['x'], data['y']]).T,\
        np.array([data[_] for _ in data_cols]).T,\
        np.array([data[_] for _ in data_errs]).T

    msk_data = dataMask(cl_data)
    ID_data, xy_data, cl_data, data_err = ID_data[msk_data],\
        xy_data[msk_data], cl_data[msk_data], data_err[msk_data]

    xy_data = MinMaxScaler().fit(xy_data).transform(xy_data)
    print("Coordinates data to normalized [0, 1]")

    # Voronoi_v4 <-- Not really needed since we are now using Voronoi on (x,y)
    # exclusively
    # cl_data, data_err = dataNorm(cl_data, data_err)

    return ID_data, xy_data, cl_data, data_err


def dataMask(data_arr, nstd=5.):
    """
    """
    msk_all = []
    # transpose is required
    for arr in data_arr.T:
        # Mask outliers (np.nan).
        med, std = np.nanmedian(arr), np.nanstd(arr)
        dmin, dmax = med - nstd * std, med + nstd * std
        msk = (arr > dmin) & (arr < dmax)
        msk_all.append(msk.data)

    msk_data = np.logical_and.reduce(msk_all)

    if (~msk_data).sum() > 0:
        print("Masked std={:.1f} outliers: {}".format(nstd, (~msk_data).sum()))

    return msk_data


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
