
import numpy as np
from . import outer
from .extras import probCnvrg
import time as t


def main(
    ID, xy, data, data_err, OL_runs, resampleFlag, PCAflag, PCAdims,
    prob_cnvrg, clust_method, otlrFlag, RK_rad, C_thresh, clust_params,
        cl_method_pars):
    """
    C_thresh : Any cluster with a smaller value will be classified as being
                composed of field stars and discarded.
    """
    start_t = t.time()

    # Set a random seed for reproducibility
    seed = np.random.randint(100000)
    print("Random seed: {}\n".format(seed))
    np.random.seed(seed)

    print("Data dimensions   : {}".format(data.shape[1]))
    print("Clustering method : {}".format(clust_method))
    print("RK rad            : {:.2f}".format(RK_rad))
    print("Threshold         : {:.1f}".format(C_thresh))

    # Initial null probabilities for all stars in the frame.
    probs_old, runs_old = np.zeros(len(ID)), 0
    probs_all = []
    for _ in range(OL_runs):
        print("\n-----------------------------------------------------------")
        print("Run {}".format(_ + 1))

        # Store all probabilities obtained in this run
        probs = outer.main(
            ID, xy, data, data_err, resampleFlag, PCAflag, PCAdims,
            clust_method, otlrFlag, RK_rad, C_thresh, clust_params,
            cl_method_pars)
        if probs:
            probs_all.append(probs)

        # DELETE
        if probs_all:
            from . import member_index
            C, P, log_MI = member_index.main(ID, np.mean(probs_all, 0))[:3]
            print("(DELETE THIS) C={:.3f}, P={:.3f}, log_MI={:.0f}".format(C, P, log_MI))

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


if __name__ == '__main__':
    main()
