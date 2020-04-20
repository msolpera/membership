
import numpy as np
from .voronoiVols import vor_2d_cummltv
from . import outer
from .extras import CCalibrate, probCnvrg
import time as t


def main(
    ID, xy, data, data_err, C_thresh=2., clust_method='agglomerative',
    method='RK', resampleFlag=False, kM_N_membs=50., kM_N_cl_max=100,
    kM_n_init=50, kM_max_iter=500, N_C_ran=100, N_runs=20, PCAflag=True,
        PCAdims='all', otlrFlag=False, perc_prob_cnvrg=.05):
    """
    C_thresh : Any cluster with a smaller value will be
                classified as being composed of field stars and discarded.
    """
    start_t = t.time()

    # Set a random seed for reproducibility
    seed = 72421 #np.random.randint(100000)
    print("Random seed: {}".format(seed))
    np.random.seed(seed)

    print("Data dimensions: {}".format(data.shape[1]))

    # Empirical CDF for the normalized 2D Voronoi areas
    vol_cummul = vor_2d_cummltv()

    if C_thresh is None:
        # Auto-calibrate the C_thresh
        from .extras import reSampleData
        from .outer import dimReduc
        data2 = reSampleData(False, data, [])
        data2 = dimReduc(data2)
        C_thresh = CCalibrate(
            xy, data2, vol_cummul, N_C_ran, kM_N_membs, kM_N_cl_max,
            kM_n_init, kM_max_iter, method)

    print("\nClustering method : {}".format(clust_method))
    print("Reject method     : {}".format(method))
    print("Threshold         : {:.3f}\n".format(C_thresh))

    # Initial null probabilities for all stars in the frame.
    probs_old, runs_old = np.zeros(len(ID)), 0
    probs_all, probs = [], []
    for _ in range(N_runs):
        print("{}.".format(_))

        # Store all probabilities obtained in this run
        probs = outer.main(
            ID, xy, data, data_err, resampleFlag, clust_method, kM_N_membs,
            kM_N_cl_max, kM_n_init, kM_max_iter, PCAflag, PCAdims, N_C_ran,
            method, vol_cummul, C_thresh, otlrFlag, probs)
        if probs:
            probs_all.append(probs)

        # DELETE
        if probs_all:
            from . import member_index
            C, P, log_MI = member_index.main(ID, np.mean(probs_all, 0))[:3]
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

    print("Minutes consumed: {:.1f}".format((t.time() - start_t) / 60.))

    return probs_mean


if __name__ == '__main__':
    main()
