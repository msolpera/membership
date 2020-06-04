
import os
from pathlib import Path
import numpy as np
from astropy.io import ascii
from astropy.stats import RipleysKEstimator
import time as t
from modules import outer
from modules.dataIO import readINI, dread, dmask, dxynorm, dwrite


def main():
    """
    """
    # Create 'output' folder if it does not exist
    Path('./output').mkdir(parents=True, exist_ok=True)

    ID_c, x_c, y_c, data_cols, data_errs, oultr_method, stdRegion_nstd,\
        rnd_seed, verbose, OL_runs, resampleFlag, PCAflag, PCAdims, GUMM_flag,\
        GUMM_perc, N_membs, clust_method, clRjctMethod, RK_rad, C_thresh,\
        cl_method_pars = readINI()

    # Process all files inside the '/input' folder
    inputfiles = readFiles()
    for file_name in inputfiles:

        if file_name.startswith("README"):
            continue
        print("\nProcessing         : {}".format(file_name))

        # Original data
        full_data, cl_ID, cl_xy, cl_data, cl_errs = dread(
            'input/' + file_name, ID_c, x_c, y_c, data_cols, data_errs,
            oultr_method, stdRegion_nstd)

        # Remove outliers
        msk_data, ID, xy, data, data_err = dmask(
            cl_ID, cl_xy, cl_data, cl_errs, oultr_method, stdRegion_nstd)

        # Normalize (x, y) data to [0, 1]
        xy01 = dxynorm(xy)

        probs_all = dataProcess(
            ID, xy01, data, data_err, rnd_seed, verbose, OL_runs, resampleFlag,
            PCAflag, PCAdims, GUMM_flag, GUMM_perc, N_membs, clust_method,
            clRjctMethod, RK_rad, C_thresh, cl_method_pars)

        if OL_runs > 1:
            # Obtain the mean of all runs. This are the final probabilities
            # assigned to each star in the frame
            probs_mean = np.mean(probs_all, 0)
        else:
            probs_mean = probs_all[0]

        # Write final data to file
        dwrite(file_name, full_data, msk_data, probs_all, probs_mean)


def dataProcess(
    ID, xy, data, data_err, rnd_seed, verbose, OL_runs, resampleFlag, PCAflag,
    PCAdims, GUMM_flag, GUMM_perc, N_membs, clust_method, clRjctMethod, RK_rad,
        C_thresh, cl_method_pars):
    """
    """
    start_t = t.time()

    # TODO this should be handled by the logging() module
    # Set print() according to the 'verbose' parameter
    if verbose == 0:
        prfl = open(os.devnull, 'w')
    else:
        prfl = None

    if clRjctMethod in ('kdetest', 'kdetestpy'):
        # Initiate empty RK_vals dictionary that will be filled by the
        # selected 'clRjctMethod'
        RK_vals, Kest = {}, None
    elif clRjctMethod == 'rkfunc':
        # Read K values from table
        RK_vals = RKDict(RK_rad)
        # Define RK test with an area of 1.
        Kest = RipleysKEstimator(area=1)

    if PCAflag:
        print("Apply PCA          : {}".format(PCAflag))
        print(" PCA N_dims        : {}".format(PCAdims))
    print("Stars per cluster  : {}".format(N_membs))
    print("Clustering method  : {}".format(clust_method))
    if cl_method_pars:
        for key, val in cl_method_pars.items():
            print(" {:<17} : {}".format(key, val))
    print("Rejection method   : {}".format(clRjctMethod))
    if clRjctMethod == 'rkfunc':
        print(" RK rad            : {:.2f}".format(RK_rad))
    print("Threshold          : {:.2f}".format(C_thresh))
    if GUMM_flag:
        print("Apply GUMM         : {}".format(GUMM_flag))
        print(" GUMM percentile   : {:.2f}".format(GUMM_perc))
    # Set a random seed for reproducibility
    if rnd_seed == 'None':
        seed = np.random.randint(100000)
    else:
        seed = int(rnd_seed)
    print("Random seed        : {}".format(seed))
    np.random.seed(seed)

    probs_all = []
    for _ in range(OL_runs):
        print("\n-----------------------------------------------------------")
        print("Run {}".format(_ + 1))

        # Store all probabilities obtained in this run
        probs, RK_vals = outer.main(
            ID, xy, data, data_err, resampleFlag, PCAflag, PCAdims, GUMM_flag,
            GUMM_perc, N_membs, clust_method, clRjctMethod, RK_rad, RK_vals,
            Kest, C_thresh, cl_method_pars, prfl)

        if probs:
            probs_all.append(probs)
        else:
            print("No probability values were obtained", file=prfl)

        p_dist = [
            (np.mean(probs_all, 0) > _).sum() for _ in (.1, .3, .5, .7, .9)]
        print("Stars with P>(.1, .3, .5, .7, .9): {}, {}, {}, {}, {}".format(
            *p_dist), file=prfl)

    elapsed = t.time() - start_t
    if elapsed > 60.:
        elapsed, ms_id = elapsed / 60., "minutes"
    else:
        ms_id = "seconds"
    print("\nTime consumed: {:.1f} {}".format(elapsed, ms_id))

    return probs_all


def RKDict(RK_rad):
    """
    Read the table with Ripley's K function pre-processed data.
    """
    # Read table with Ripley's K stored data
    RK_data = ascii.read("modules/K_table.dat")

    # Generate the final dictionary with 'N' as keys, and (mean, std) as values
    keys = RK_data['N']
    vals = np.array(
        [RK_data['mean_' + str(RK_rad)].data,
         RK_data['std_' + str(RK_rad)].data]).T
    RK_vals = dict(zip(keys, vals))

    return RK_vals


def readFiles():
    """
    Read files from the input folder
    """
    return [arch.name for arch in Path('input').iterdir() if arch.is_file()]


if __name__ == '__main__':
    main()
