
import numpy as np
from .voronoiVols import vor_2d_cummltv
from .extras import kMeansClust, rjctField, outlierRjct, CCalibrate,\
    probCnvrg, reSampleData
from . import member_index


def main(
    ID, xy, data, data_err, C_thresh=2, method='AD', resampleFlag=False,
    kM_N_membs=50., kM_N_cl_max=100, kM_n_init=20, kM_max_iter=100, N_C_ran=100,
        N_runs=10, perc_prob_cnvrg=.05):
    """
    C_thresh : maximum AD|KS value. Any cluster with a smaller value will be
                classified as being composed of field stars and discarded.

    For a two-sample test, the AD (critical) values:
       0.325,  1.226,  1.961,  2.718,  3.752, 4.592, 6.546
    correspond to significance levels of (in %):
          25,     10,      5,    2.5,      1,   0.5,   0.1

    A larger threshold improves the Purity (P) at the expense of the
    Completeness (C), and vice-versa.
    """

    # Set a random seed for reproducibility
    seed = 72421 #np.random.randint(100000)
    print("Random seed: {}".format(seed))
    np.random.seed(seed)

    print("Data dimensions: {}".format(data.shape[1]))

    # Empirical CDF for the normalized 2D Voronoi areas
    vol_cummul = vor_2d_cummltv()

    # if C_thresh is None:
    #     # Auto-calibrate the C_thresh
    #     C_thresh = CCalibrate(
    #         xy, data, vol_cummul, N_C_ran, kM_N_membs, kM_N_cl_max, kM_n_init,
    #         kM_max_iter, method)

    print("\nMethod: {}, Threshold: {:.3f}\n".format(method, C_thresh))

    # Initial null probabilities for all stars in the frame.
    probs_old, runs_old = np.zeros(len(ID)), 0
    probs_all = []
    for C_thresh in np.linspace(0., 3., 20): # range(N_runs):
        print("{}.".format(C_thresh))

        # Make a copy of the original data to avoid over-writing it
        clust_ID, clust_xy = np.array(list(ID)), np.array(list(xy))

        # Re-sample the data using its uncertainties?
        clust_data = reSampleData(resampleFlag, data, data_err, True)

        # Keep going until all the "fake clusters" are rejected
        while True:

            # Obtain all the clusters in the input data using kMeans
            clusts_msk = kMeansClust(
                clust_data, kM_N_membs, kM_N_cl_max, kM_n_init,
                kM_max_iter)

            # Masks of rejected "field" clusters
            C_masks = rjctField(
                N_C_ran, method, vol_cummul, C_thresh, clust_xy, clusts_msk)

            if not C_masks:
                # No more clusters to reject
                break

            # Combine all the masks
            msk_all = np.logical_and.reduce(C_masks)

            # Applying this mask leaves too few stars --> Break
            if clust_data[msk_all].shape[0] < 10:
                print(" N_stars<10 after rejecting clusters. Breaking")
                break

            # Remove stars identified as field stars from the frame and move on
            # to the next iteration
            clust_ID, clust_xy, clust_data = clust_ID[msk_all],\
                clust_xy[msk_all], clust_data[msk_all]

        # Mark all the stars that survived in 'clust_ID' as members assigning
        # a probability of '1'. All others are field stars and are assigned
        # a probability of '0'.
        probs = []
        for st in ID:
            if st in clust_ID:
                probs.append(1.)
            else:
                probs.append(0.)

        probs = outlierRjct(clust_xy, probs, False)

        # Store all probabilities obtained in this run
        probs_all.append(probs)

        C, P, log_MI = member_index.main(ID, probs)[:3]
        print("C={:.3f}, P={:.3f}, log_MI={:.0f}".format(C, P, log_MI))

        # Break out if/when probabilities converge
        cnvrg_flag, probs_old, runs_old = probCnvrg(
            probs_all, probs_old, perc_prob_cnvrg, runs_old)
        if cnvrg_flag:
            print("Probabilities converged to {}%. Breaking".format(
                perc_prob_cnvrg))
            break

    # Obtain the mean of all runs. This is the final probabilities assigned
    # to each star in the frame
    probs_mean = np.mean(probs_all, 0)

    return probs_mean


if __name__ == '__main__':
    main()
