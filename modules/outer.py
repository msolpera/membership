
import numpy as np
from . import inner
from .extras import reSampleData, outlierRjct
from sklearn.decomposition import PCA


def main(
    ID, xy, data, data_err, resampleFlag, PCAflag, PCAdims, clust_method,
        otlrFlag, C_thresh, unif_method, RK_rad, clust_params, vol_cummul):
    """
    """

    # Make a copy of the original data to avoid over-writing it
    clust_ID, clust_xy = np.array(list(ID)), np.array(list(xy))

    # Re-sample the data using its uncertainties?
    clust_data = reSampleData(resampleFlag, data, data_err)

    # Apply PCA and features reduction
    clust_data = dimReduc(clust_data, PCAflag, PCAdims)

    nostars_flag = False
    # Keep calling the inner loop until all the "fake clusters" are rejected
    while True:

        C_masks, N_clusts = inner.main(
            clust_xy, clust_data, clust_method, otlrFlag, C_thresh,
            unif_method, RK_rad, clust_params, vol_cummul)

        if not C_masks:
            # No more clusters to reject
            break

        # Combine all the masks
        msk_all = np.logical_and.reduce(C_masks)

        # This mask leaves too few stars and too many clusters --> Break
        if clust_data[msk_all].shape[0] < int(clust_params['N_membs']):
            # If there are more than 4 (HARDCODED) clusters defined at this
            # point, reject this run.
            if N_clusts > 4:
                nostars_flag = True
            print(" N stars<{:.0f} Breaking".format(int(
                clust_params['N_membs'])))
            break

        # Remove stars identified as field stars from the frame and move on
        # to the next iteration
        clust_ID, clust_xy, clust_data = clust_ID[msk_all],\
            clust_xy[msk_all], clust_data[msk_all]

    probs = []
    if nostars_flag is False:
        # Mark all the stars that survived in 'clust_ID' as members assigning
        # a probability of '1'. All others are field stars and are assigned
        # a probability of '0'.
        for st in ID:
            if st in clust_ID:
                probs.append(1.)
            else:
                probs.append(0.)

        # DELETE
        from . import member_index
        C, P, log_MI = member_index.main(ID, probs)[:3]
        print("C={:.3f}, P={:.3f}, log_MI={:.0f}".format(C, P, log_MI))

        probs = outlierRjct(clust_xy, probs, otlrFlag)

    return probs


def dimReduc(cl_data, PCAflag, PCAdims):
    """
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
