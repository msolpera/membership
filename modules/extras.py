
import numpy as np
from astropy.stats import sigma_clipped_stats
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import LocalOutlierFactor
from sklearn.ensemble import IsolationForest
from .voronoiVols import voronoi_volumes
from . import compFunc


def CCalibrate(
    clust_xy, cl_data, clust_method, vol_cummul, N_C_ran, kM_N_membs,
        kM_N_cl_max, kM_n_init, kM_max_iter, method):
    """
    """

    # vol_d = voronoi_volumes(clust_xy)
    # import matplotlib.pyplot as plt
    # from .voronoiVols import vor_2d_poisson
    # xx =  np.linspace(0., 5., 1000)
    # yy = vor_2d_poisson(xx)
    # plt.hist(vol_d/vol_d.mean(), 100, alpha=.25, color='b',density=True)
    # plt.plot(xx,yy)

    print("Calibrating C threshold")
    from .inner import clustAlgor

    # Obtain all the clusters in the input data using kMeans
    clusts_msk = clustAlgor(
        cl_data, clust_method, kM_N_membs, kM_N_cl_max, kM_n_init, kM_max_iter,
        [])

    C_vals = []
    for i, cl_msk in enumerate(clusts_msk):
        if cl_msk.sum() >= 5:
            xy = clust_xy[cl_msk]
            vol_d = voronoi_volumes(xy)
            C_vals.append(compFunc.main(
                method, N_C_ran, xy, vol_cummul, vol_d))

    # Using this median gives good results.
    mean, median, std = sigma_clipped_stats(C_vals)
    print(mean, median, std)
    # print(np.median(C_vals), np.mean(C_vals), np.std(C_vals))
    C_thresh = np.mean(C_vals)
    import pdb; pdb.set_trace()  # breakpoint 0c2ef80c //


    return C_thresh


def reSampleData(resampleFlag, data, data_err, standard_scale=True):
    """
    Re-sample the data given its uncertainties using a normal distribution
    """
    if resampleFlag:
        # Gaussian random sample
        grs = np.random.normal(0., 1., data.shape[0])
        sampled_data = data + grs[:, np.newaxis] * data_err
    else:
        sampled_data = np.array(list(data))

    if standard_scale:
        sampled_data = StandardScaler().fit(sampled_data).transform(
            sampled_data)

    return sampled_data


def outlierRjct(clust_xy, probs, outrjctFlag=True, n_neighbors=50):
    """
    Reject outliers in the (x, y) coordinates space using the 'Unsupervised
    Outlier Detection using Local Outlier Factor (LOF)' method.
    """

    if outrjctFlag is False:
        print("No outlier re-classification")

    elif clust_xy.shape[0] < n_neighbors:
        print("Not enough stars (<{}) to perform outlier rejection".format(
            n_neighbors))

    else:
        # Predict outliers
        y_pred = LocalOutlierFactor(n_neighbors=n_neighbors).fit_predict(
            clust_xy)
        # y_pred = IsolationForest().fit_predict(clust_xy)

        j, Nr = 0, 0
        for i, p in enumerate(probs):
            # This star was classified as member
            if p > .5:
                # But rejected as outlier
                if y_pred[j] == -1:
                    # Classify as field
                    probs[i] = 0.
                    Nr += 1
                j += 1

        print("Stars re-classified as field: {}".format(Nr))

    return probs


def probCnvrg(probs_all, probs_old, perc_prob_cnvrg, runs_old, min_runs=3):
    """
    Check if probabilities converged to within 'perc_prob_cnvrg'%.

    Break only if a minimum of 'min_runs' consecutive runs have been processed.
    """
    cnvrg_flag = False

    if not probs_all:
        return cnvrg_flag, probs_old, runs_old

    # Average all probabilities.
    prob_avrg = np.mean(probs_all, 0)

    # if runs_old > min_runs:
    #     for _ in np.linspace(0.01, .5, 50):
    #         if np.allclose(probs_old, prob_avrg, _):
    #             print("Relative tolerance for probabilities: {:.0f}%".format(
    #                 100. * _))
    #             break

    # Check if probabilities changed less than 'perc_prob_cnvrg'% with
    # respect to the previous iteration.
    if np.allclose(probs_old, prob_avrg, perc_prob_cnvrg):
        if runs_old > min_runs:
            # Arrays are equal within tolerance
            cnvrg_flag = True
        runs_old += 1
    else:
        # Store new array in old one and proceed to new iteration.
        probs_old = prob_avrg
        # Reset
        runs_old = 0

    return cnvrg_flag, probs_old, runs_old
