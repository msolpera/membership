
import numpy as np
from scipy.stats import gaussian_kde
from . import clustAlgor


def main(
    clust_xy, clust_data, N_membs, clust_method, clRjctMethod, RK_rad,
        RK_vals, Kest, C_thresh, cl_method_pars, prfl, N_cl_max=1000):
    """
    Perform the inner loop: cluster --> reject
    """

    print("  Performing clustering on array of shape ({}, {})".format(
        *clust_data.shape), file=prfl)

    # Number of clusters: min is 2, max is N_cl_max
    n_clusters = max(2, min(int(clust_data.shape[0] / N_membs), N_cl_max))

    # Obtain all the clusters in the input data using kMeans
    if clust_method == 'Voronoi':
        labels = clustAlgor.voronoi(clust_data, n_clusters)
    elif clust_method == 'rkde':
        labels = clustAlgor.RKDE(clust_data, n_clusters)
    # Not fully implemented yet
    elif clust_method[:4] != 'pycl':
        labels = clustAlgor.sklearnMethods(
            clust_method, cl_method_pars, clust_data, n_clusters)
    else:
        labels = clustAlgor.pycl(clust_data, n_clusters)

    # Separate the labels that point to each cluster found
    clusts_msk = [(labels == _) for _ in range(labels.min(), labels.max() + 1)]

    # import matplotlib.pyplot as plt
    # for cl in clusts_msk:
    #     plt.scatter(*clust_data[cl].T, alpha=.5)

    print("  Identified {} clusters".format(len(clusts_msk)), file=prfl)

    # Reject "field" clusters and return their stored masks
    C_masks, RK_vals = rjctField(
        clust_xy, clusts_msk, clRjctMethod, RK_rad, RK_vals, Kest, C_thresh,
        prfl)

    return C_masks, RK_vals, len(clusts_msk)


def rjctField(
    clust_xy, clusts_msk, clRjctMethod, RK_rad, RK_vals, Kest, C_thresh, prfl,
        minStars=5):
    """
    Iterate over all the clusters defined by the clustering algorithm, storing
    those masks that point to clusters that survived the test with uniform
    random 2D fields in (x, y).

    C_thresh : Any cluster with a smaller value will be classified as being
                composed of field stars and discarded.
    minStars: int
        If the cluster contains less than this many stars, skip the test and
        mark its stars as field stars
    """

    C_masks = []
    # For each cluster found by the clustering method, check if it is composed
    # of field stars or actual cluster members, using their (x, y)
    # distribution, and Ripley's K test.
    for i, cl_msk in enumerate(clusts_msk):

        if cl_msk.sum() < minStars:
            C_s = -np.inf
        else:
            C_s, RK_vals = rejctVal(
                clust_xy[cl_msk], clRjctMethod, RK_rad, RK_vals, Kest)

        # Smaller C_s values point to samples that come from a uniform
        # distribution in (x, y), i.e., a "cluster" made up of field
        # stars. Hence, we keep as "true" clusters those with C_s values
        # larger than C_thresh.
        if C_s >= C_thresh:
            # Store mask that points to stars that should be *kept*
            print("   Cluster {} survived (C_s={:.2f}), N={}".format(
                i, C_s, cl_msk.sum()), file=prfl)
            C_masks.append(cl_msk)

    return C_masks, RK_vals


def rejctVal(xy, clRjctMethod, RK_rad, RK_vals, Kest):
    """
    Test how similar this cluster's (x, y) distribution is compared
    to a uniform random distribution.
    """

    if clRjctMethod == 'rkfunc':
        # Test how similar this cluster's (x, y) distribution is compared
        # to a uniform random distribution using Ripley's K.
        # https://stats.stackexchange.com/a/122816/10416

        # Read stored value from table.
        try:
            mean, std = RK_vals[xy.shape[0]]
        except KeyError:
            # If the number of stars in this cluster is larger than the largest
            # one stored, used the values for this largest value.
            mean, std = list(RK_vals.values())[-1]

        C_s = (Kest(xy, (RK_rad,))[0] - mean) / std

    elif clRjctMethod == "kdetest":

        from rpy2.robjects import r
        from rpy2.robjects import numpy2ri
        from rpy2.robjects.packages import importr
        importr('MASS')
        numpy2ri.activate()

        N = xy.shape[0]
        xmin, ymin = xy.min(0)
        xmax, ymax = xy.max(0)
        xrng, yrng = xmax - xmin, ymax - ymin

        r.assign('minX', xmin)
        r.assign('maxX', xmax)
        r.assign('minY', ymin)
        r.assign('maxY', ymax)
        r.assign('xrng', xrng)
        r.assign('yrng', yrng)
        r.assign('nstars', N)

        # Read stored value from table.
        try:
            mean, std = RK_vals[xy.shape[0]]
        except KeyError:
            dist_u = []
            for _ in range(100):
                r('dataX <- runif(nstars, 0, xrng)')
                r('dataY <- runif(nstars, 0, yrng)')
                r('kde2dmap <- kde2d(dataX, dataY, n=50, lims=c(0, xrng, 0, yrng))')
                kde2dmap = np.array(list(r('kde2dmap'))[2]).flatten()
                dist_u.append(
                    (kde2dmap.max() - kde2dmap.mean()) / kde2dmap.std())

            mean, std = np.mean(dist_u), np.std(dist_u)
            RK_vals[N] = mean, std

        rx = r.matrix(xy.T[0])
        ry = r.matrix(xy.T[1])
        r.assign('dataX', rx)
        r.assign('dataY', ry)
        r('kde2dmap <- kde2d(dataX, dataY, n=50, lims=c(minX, maxX, minY, maxY))')
        kde2dmap = np.array(list(r('kde2dmap'))[2]).flatten()
        dist_d = (kde2dmap.max() - kde2dmap.mean()) / kde2dmap.std()

        C_s = (dist_d - mean) / std

    elif clRjctMethod == "kdetestpy":
        N = xy.shape[0]
        xmin, ymin = 0., 0.
        xmax, ymax = xy.max(0)
        xx, yy = np.mgrid[xmin:xmax:50j, ymin:ymax:50j]
        positions = np.vstack([xx.ravel(), yy.ravel()])

        # Read stored value from table.
        try:
            mean, std = RK_vals[xy.shape[0]]
        except KeyError:
            dist_u = []
            for _ in range(100):
                # Generate random uniform 2D distribution
                xy_u = np.random.uniform((xmin, ymin), (xmax, ymax), (N, 2))
                kde = gaussian_kde(xy_u.T)
                kde2dmap = kde.evaluate(positions)
                dist_u.append(
                    (kde2dmap.max() - kde2dmap.mean()) / kde2dmap.std())

            mean, std = np.mean(dist_u), np.std(dist_u)
            RK_vals[N] = mean, std

        # Evaluate subset
        # from scipy.stats import iqr
        # bw = 1.06 * np.min([xy.std(None), iqr(xy, None) / 1.34], None) *\
        #     N**(-1. / 5.)
        # kde = gaussian_kde(xy.T, bw_method=bw / xy.T.std(ddof=1))
        xmin, ymin = xy.min(0)
        xmax, ymax = xy.max(0)
        xx, yy = np.mgrid[xmin:xmax:50j, ymin:ymax:50j]
        positions = np.vstack([xx.ravel(), yy.ravel()])
        kde = gaussian_kde(xy.T)
        kde2dmap = kde.evaluate(positions)
        dist_d = (kde2dmap.max() - kde2dmap.mean()) / kde2dmap.std()

        C_s = (dist_d - mean) / std

    return C_s, RK_vals


# C_vals = []
# if method == '1':
#     # This methods uses always the original xy distribution and compares it
#     # with a uniform distribution
#     K_d = Kest(xy, (RK_rad,), mode=RK_mode)
#     for _ in range(N_C_ran):
#         xy_u = np.random.uniform(0., 1., xy.shape)
#         K_u = Kest(xy_u, (RK_rad,), mode=RK_mode)
#         C_vals.append(K_d / K_u)
#     # Percentage of times that K_d >= K_u (i.e., the xy array was
#     # identified as a cluster)
#     C_s = (np.array(C_vals) >= 1.).sum() / N_C_ran
# elif method == '2':
#     # This methods bootstraps the original xy distribution and compares it
#     # with a uniform distribution
#     for _ in range(N_C_ran):
#         xy_u = np.random.uniform(0., 1., xy.shape)
#         K_u = Kest(xy_u, (RK_rad,), mode=RK_mode)
#         # Bootstrap original xy distribution
#         msk = np.random.choice(xy.shape[0], xy.shape[0])
#         xy_b = xy[msk]
#         K_d = Kest(xy_b, (RK_rad,), mode=RK_mode)
#         C_vals.append(K_d / K_u)
#     # Percentage of times that K_d >= K_u (i.e., the xy array was
#     # identified as a cluster)
#     C_s = (np.array(C_vals) >= 1.).sum() / N_C_ran
