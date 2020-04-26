
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.neighbors import LocalOutlierFactor
from sklearn.ensemble import IsolationForest


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


def dimReduc(cl_data, PCAflag, PCAdims):
    """
    Perform PCA and feature reduction

    all: use all available dimensions
    """
    if PCAflag:
        if PCAdims == 'all':
            PCAdims = cl_data.shape[1]
        else:
            PCAdims = int(PCAdims)

        pca = PCA(n_components=PCAdims)
        cl_data_pca = pca.fit(cl_data).transform(cl_data)
        print("Selected N={} PCA features".format(PCAdims))
    else:
        cl_data_pca = cl_data

    return cl_data_pca


def outlierRjct(clust_xy, probs, outrjctFlag=True, n_neighbors=50):
    """
    Reject outliers in the (x, y) coordinates space using the 'Unsupervised
    Outlier Detection using Local Outlier Factor (LOF)' method.
    """

    if outrjctFlag is False:
        # print("No outlier re-classification")
        pass

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
