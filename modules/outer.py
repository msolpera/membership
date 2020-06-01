
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from . import inner
from .GUMM import GUMMtrain


def main(
    ID, xy, data, data_err, resampleFlag, PCAflag, PCAdims, GUMM_flag,
    GUMM_perc, N_membs, clust_method, clRjctMethod, RK_rad, RK_vals,
        Kest, C_thresh, cl_method_pars, prfl):
    """
    Perform the outer loop: inner loop until all "fake" clusters are rejected
    """

    # Make a copy of the original data to avoid over-writing it
    clust_ID, clust_xy = np.array(list(ID)), np.array(list(xy))

    # Re-sample the data using its uncertainties?
    clust_data = reSampleData(resampleFlag, data, data_err)
    # Apply PCA and features reduction
    clust_data = dimReduc(clust_data, PCAflag, PCAdims)
    # clust_data = np.array(list(data))

    # Keep calling the inner loop until all the "fake clusters" are rejected
    _iter = 1
    while True:
        print("\n Iteration {}".format(_iter), file=prfl)
        _iter += 1

        C_masks, RK_vals, N_clusts = inner.main(
            clust_xy, clust_data, N_membs, clust_method, clRjctMethod, RK_rad,
            RK_vals, Kest, C_thresh, cl_method_pars, prfl)

        # No clusters were rejected in this iteration. Break
        if N_clusts == len(C_masks):
            print(" All clusters survived, N={}".format(
                clust_xy.shape[0]), file=prfl)
            break

        # Combine all the masks using a logical OR
        msk_all = np.logical_or.reduce(C_masks)

        # Applying 'msk_all' results in too few stars distributed in too many
        # clusters. Break
        if clust_data[msk_all].shape[0] < int(N_membs):
            print(" N_stars<{:.0f} Breaking".format(int(N_membs)), file=prfl)
            break

        # Keep only stars identified as members and move on to the next
        # iteration
        clust_ID, clust_xy, clust_data = clust_ID[msk_all],\
            clust_xy[msk_all], clust_data[msk_all]
        print(" A total of {} stars survived in {} clusters".format(
            msk_all.sum(), len(C_masks)), file=prfl)

    probs = []
    # Mark all the stars that survived in 'clust_ID' as members assigning
    # a probability of '1'. All others are field stars and are assigned
    # a probability of '0'.
    for i, st in enumerate(ID):
        if st in clust_ID:
            probs.append(1.)
        else:
            probs.append(0.)

    if GUMM_flag:
        # Perform a cleaning on the final list of stars selected as members.
        probs = GUMMtrain(GUMM_perc, clust_xy, probs)

    return probs, RK_vals


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
        print("Standard scale: removed mean and scaled to unit variance")
        sampled_data = StandardScaler().fit(sampled_data).transform(
            sampled_data)

    return sampled_data


def dimReduc(cl_data, PCAflag, PCAdims):
    """
    Perform PCA and feature reduction
    """
    if PCAflag:
        pca = PCA(n_components=PCAdims)
        cl_data_pca = pca.fit(cl_data).transform(cl_data)
        print("Selected N={} PCA features".format(PCAdims))
    else:
        cl_data_pca = cl_data

    return cl_data_pca
