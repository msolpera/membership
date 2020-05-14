
import os
import numpy as np
from astropy.io import ascii
from astropy.stats import RipleysKEstimator
from . import outer
from .extras import probCnvrg
import time as t


def main(
    ID, xy, data, data_err, verbose, OL_runs, resampleFlag, PCAflag, PCAdims,
    prob_cnvrg, clust_method, RK_rad, RK_mode, C_thresh, clust_params,
        cl_method_pars):
    """
    C_thresh : Any cluster with a smaller value will be classified as being
                composed of field stars and discarded.
    """
    start_t = t.time()

    # Set print() according to the 'verbose' parameter
    if verbose == 0:
        prfl = open(os.devnull, 'w')
    else:
        prfl = None

    RK_vals = RKDict(RK_rad)

    # TODO: make this more general?
    # Define the (x, y) area with sides [0, 1]
    Kest = RipleysKEstimator(area=1)  # , x_max=1, y_max=1, x_min=0, y_min=0)

    # Set a random seed for reproducibility
    seed = np.random.randint(100000)
    print("Random seed: {}\n".format(seed))
    np.random.seed(seed)

    print("Data dimensions   : {}".format(data.shape[1]))
    print("Stars per cluster : {}".format(clust_params['N_membs']))
    print("Clustering method : {}".format(clust_method))
    if cl_method_pars:
        for key, val in cl_method_pars.items():
            print(" {:<16} : {}".format(key, val))
    print("RK rad            : {:.2f}".format(RK_rad))
    # print("RK mode           : {}".format(RK_mode))
    print("Threshold         : {:.2f}".format(C_thresh))

    clust_ID = np.array(list(ID))

    # Initial null probabilities for all stars in the frame.
    probs_old, runs_old = np.zeros(len(ID)), 0
    probs_all = []
    for _ in range(OL_runs):
        print("\n-----------------------------------------------------------")
        print("Run {}".format(_ + 1))

        # Store all probabilities obtained in this run
        probs = outer.main(
            ID, xy, data, data_err, resampleFlag, PCAflag, PCAdims,
            clust_method, RK_rad, RK_vals, Kest, C_thresh,
            clust_params, cl_method_pars, prfl)

        if probs:
            probs_all.append(probs)

        # DELETE
        if probs_all:
            from . import member_index
            C, P, log_MI = member_index.main(
                clust_ID, np.mean(probs_all, 0))[:3]
            print("C={:.3f}, P={:.3f}, log_MI={:.0f}".format(
                C, P, log_MI))
        # DELETE

        p_dist = [
            (np.mean(probs_all, 0) > _).sum() for _ in (.1, .3, .5, .7, .9)]
        print("Stars with P>(.1, .3, .5, .7, .9): {}, {}, {}, {}, {}".format(
            *p_dist))

        # Break out if/when probabilities converge
        cnvrg_flag, probs_old, runs_old = probCnvrg(
            probs_all, probs_old, prob_cnvrg, runs_old)
        if cnvrg_flag:
            print("Probabilities converged to {}%. Breaking".format(
                prob_cnvrg))
            break

    # Obtain the mean of all runs. This is the final probabilities assigned
    # to each star in the frame
    probs_mean = np.mean(probs_all, 0)

    print("Minutes consumed: {:.1f}".format((t.time() - start_t) / 60.))

    return probs_mean


def RKDict(RK_rad):
    """
    Read the table with Ripley's K function pre-processed data.
    """
    # Read table with Ripley's K stored data
    RK_data = ascii.read("modules/RK_data.dat")

    # Generate the final dictionary with 'N' as keys, and (mean, std) as values
    keys = RK_data['N']
    vals = np.array(
        [RK_data[str(RK_rad) + '_m'].data, RK_data[str(RK_rad) + '_s'].data]).T
    RK_vals = dict(zip(keys, vals))

    return RK_vals


if __name__ == '__main__':
    main()
