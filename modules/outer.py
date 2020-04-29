
import numpy as np
from . import inner
from .extras import reSampleData, dimReduc, outlierRjct


def main(
    ID, xy, data, data_err, resampleFlag, PCAflag, PCAdims, clust_method,
    otlrFlag, RK_rad, RK_mode, C_thresh, clust_params, cl_method_pars,
        probs_outer):
    """
    Perform the outer loop: inner loop until all "fake" clusters are rejected
    """

    # Make a copy of the original data to avoid over-writing it
    clust_ID, clust_xy = np.array(list(ID)), np.array(list(xy))

    # Re-sample the data using its uncertainties?
    clust_data = reSampleData(resampleFlag, data, data_err)

    # Apply PCA and features reduction
    clust_data = dimReduc(clust_data, PCAflag, PCAdims)

    # Keep calling the inner loop until all the "fake clusters" are rejected
    _iter, nostars_flag = 1, False
    while True:
        print("\n Iteration {}".format(_iter))
        _iter += 1

        C_masks, N_clusts = inner.main(
            clust_xy, clust_data, clust_method, RK_rad, RK_mode, C_thresh,
            clust_params, cl_method_pars, probs_outer)

        # No clusters were rejected in this iteration. Break
        if N_clusts == len(C_masks):
            print(" All clusters ({} stars) survived".format(
                clust_xy.shape[0]))
            break

        # Combine all the masks using a logical OR
        msk_all = np.logical_or.reduce(C_masks)

        # Applying 'msk_all' results in too few stars distributed in too many
        # clusters. Break
        if clust_data[msk_all].shape[0] < int(clust_params['N_membs']):
            # If there are more than 4 (HARDCODED) clusters defined at this
            # point, reject this run.
            if N_clusts > 4:
                nostars_flag = True
            print(" N_stars<{:.0f} Breaking".format(int(
                clust_params['N_membs'])))
            break

        # Keep only stars identified as members and move on to the next
        # iteration
        clust_ID, clust_xy, clust_data = clust_ID[msk_all],\
            clust_xy[msk_all], clust_data[msk_all]
        print(" A total of {} stars survived in {} clusters".format(
            msk_all.sum(), len(C_masks)))

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
        print("\n(DELETE THIS) C={:.3f}, P={:.3f}, log_MI={:.0f}".format(C, P, log_MI))

        probs = outlierRjct(clust_xy, probs, otlrFlag)

    return probs
