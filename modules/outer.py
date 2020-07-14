
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from . import inner
from .GUMM import GUMMProbs
from .GUMMExtras import GUMMProbCut, lowCIGUMMClean


def loop(
    ID, xy, data, data_err, resampleFlag, PCAflag, PCAdims, GUMM_flag,
    GUMM_perc, N_membs, clust_method, clRjctMethod, Kest, C_thresh,
        cl_method_pars, prfl, prob_GUMM, KDE_vals):
    """
    Perform the outer loop: inner loop until all "fake" clusters are rejected
    """

    # Make a copy of the original data to avoid over-writing it
    clust_ID, clust_xy = np.array(list(ID)), np.array(list(xy))

    # Re-sample the data using its uncertainties?
    clust_data = reSampleData(resampleFlag, data, data_err, prfl)
    # Apply PCA and features reduction
    clust_data = dimReduc(clust_data, PCAflag, PCAdims, prfl)

    # Keep calling the inner loop until all the "fake clusters" are rejected
    _iter = 1
    while True:
        print("\n IL iteration {}".format(_iter), file=prfl)
        _iter += 1

        C_masks, KDE_vals, N_clusts = inner.main(
            clust_xy, clust_data, N_membs, clust_method, clRjctMethod,
            KDE_vals, Kest, C_thresh, cl_method_pars, prfl)

        # No clusters were rejected in this iteration. Break
        if N_clusts == len(C_masks):
            print(" All clusters survived, N={}".format(
                clust_xy.shape[0]), file=prfl)
            break

        # Combine all the masks using a logical OR
        msk_all = np.logical_or.reduce(C_masks)

        # Applying 'msk_all' results in too few stars. Break
        if clust_data[msk_all].shape[0] < N_membs:
            print(" N_stars<{:.0f} Breaking".format(N_membs), file=prfl)
            break

        # Keep only stars identified as members and move on to the next
        # iteration
        clust_ID, clust_xy, clust_data = clust_ID[msk_all],\
            clust_xy[msk_all], clust_data[msk_all]
        print(" A total of {} stars survived in {} clusters".format(
            msk_all.sum(), len(C_masks)), file=prfl)

        # Clean using GUMM
        if GUMM_flag:
            gumm_p = GUMMProbs(clust_xy, prfl)
            prob_cut = GUMMProbCut(GUMM_perc, gumm_p, prob_GUMM)
            # Mark all stars as members
            probs_cl = np.ones(len(clust_xy))
            # Mark as non-members those below 'prob_cut'
            probs_cl[gumm_p <= prob_cut] = 0.

            # Keep only member stars for the next run (if enough stars remain
            # in the list)
            msk = probs_cl > 0.
            if msk.sum() > N_membs:
                clust_ID, clust_xy, clust_data = clust_ID[msk], clust_xy[msk],\
                    clust_data[msk]
                print("GUMM analysis: reject {} stars as non-members".format(
                    msk.sum()), file=prfl)

    # Mark all the stars that survived in 'clust_ID' as members assigning
    # a probability of '1'. All others are field stars and are assigned
    # a probability of '0'.
    cl_probs = np.zeros(len(ID))
    for i, st in enumerate(ID):
        if st in clust_ID:
            cl_probs[i] = 1.

    # Perform a final cleaning on the list of stars selected as members.
    # Use the last list of coordinates and IDs from the inner loop.
    if GUMM_flag:
        # This is only ever used for *very* low contaminated clusters.
        cl_probs = lowCIGUMMClean(
            N_membs, GUMM_perc, ID, cl_probs, clust_ID, clust_xy, prfl, prob_GUMM)

    return list(cl_probs), KDE_vals


def reSampleData(resampleFlag, data, data_err, prfl, standard_scale=True):
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
        print(
            "Standard scale: removed mean and scaled to unit variance",
            file=prfl)
        sampled_data = StandardScaler().fit(sampled_data).transform(
            sampled_data)

    return sampled_data


def dimReduc(cl_data, PCAflag, PCAdims, prfl):
    """
    Perform PCA and feature reduction
    """
    if PCAflag:
        pca = PCA(n_components=PCAdims)
        cl_data_pca = pca.fit(cl_data).transform(cl_data)
        print("Selected N={} PCA features".format(PCAdims), file=prfl)
        var_r = ["{:.2f}".format(_) for _ in pca.explained_variance_ratio_]
        print(" Variance ratio: ", ", ".join(var_r), file=prfl)
    else:
        cl_data_pca = cl_data

    return cl_data_pca
