
import warnings
import numpy as np
from scipy.stats import gaussian_kde
from . import clustAlgor


def main(
    clust_xy, clust_data, N_membs, clust_method, clRjctMethod,
    KDE_vals, Kest, C_thresh, cl_method_pars, prfl, N_cl_max=10000,
        N_st_max=20000):
    """
    Perform the inner loop: cluster --> reject

    N_cl_max : maximum number of clusters allowed
    N_st_max : maximum numbers of stars for the Voronoi and kNNdens algorithms.
               Above this threshold, the 'cdist' function is avoided.
    """

    print("  Performing clustering on array of shape ({}, {})".format(
        *clust_data.shape), file=prfl)

    # Number of clusters: min is 2, max is N_cl_max
    n_clusters = max(2, int(clust_data.shape[0] / N_membs))
    if n_clusters > N_cl_max:
        print("  Too many clusters. Capping at {}".format(N_cl_max), file=prfl)
        n_clusters = N_cl_max

    if clust_method == 'Voronoi':
        labels = clustAlgor.voronoi(clust_data, N_membs, n_clusters, N_st_max)

    # Obtain all the clusters in the input data
    elif clust_method == 'kNNdens':
        labels = clustAlgor.kNNdens(
            clust_data, cl_method_pars, N_membs, n_clusters, N_st_max)

    elif clust_method == 'rkmeans':
        labels = clustAlgor.RKmeans(clust_data, n_clusters)

    # scikit-learn methods
    elif clust_method[:4] != 'pycl':
        labels = clustAlgor.sklearnMethods(
            clust_method, cl_method_pars, clust_data, n_clusters)

    # TODO Not fully implemented yet
    else:
        labels = clustAlgor.pycl(clust_data, n_clusters)

    # import matplotlib.pyplot as plt
    # for cl in clusts_msk:
    #     plt.scatter(*clust_data[cl].T[:2], alpha=.25)

    N_clusts = len(set(labels))
    print("  Identified {} clusters".format(N_clusts), file=prfl)

    # Reject "field" clusters and return their stored masks
    C_masks, N_survived, KDE_vals = rjctField(
        clust_xy, labels, clRjctMethod, KDE_vals, Kest, C_thresh,
        prfl)

    return N_clusts, C_masks, N_survived, KDE_vals


def rjctField(
    clust_xy, labels, clRjctMethod, KDE_vals, Kest, C_thresh, prfl,
        minStars=5):
    """
    Iterate over all the clusters defined by the clustering algorithm, storing
    those masks that point to clusters that survived the test with uniform
    random 2D fields in (x, y).

    C_thresh : Any cluster with a smaller value will be classified as being
                composed of field stars and discarded.
    minStars: int
        If the cluster contains less than this many stars, skip the test and
        mark its stars as field stars.
    """

    C_masks, N_survived = [], 0
    # For each cluster found by the clustering method, check if it is composed
    # of field stars or actual cluster members, using their (x, y)
    # distribution, and Ripley's K test.
    # for i, cl_msk in enumerate(clusts_msk):
    for i in range(labels.min(), labels.max() + 1):
        # Separate stars assigned to this label
        cl_msk = labels == i

        if cl_msk.sum() < minStars:
            C_s = -np.inf
        else:
            C_s, KDE_vals = rejctVal(
                clust_xy[cl_msk], clRjctMethod, KDE_vals, Kest)

        # Smaller C_s values point to samples that come from a uniform
        # distribution in (x, y), i.e., a "cluster" made up of field
        # stars. Hence, we keep as "true" clusters those with C_s values
        # larger than C_thresh.

        if clRjctMethod == 'rkfunc':
            # 1% critical value. From Dixon (2001), 'Ripley's K function'
            C_thresh = 1.68 / cl_msk.sum()

        if C_s >= C_thresh:
            # Store mask that points to stars that should be *kept*
            print("   Cluster {} survived (C_s={:.2f}), N={}".format(
                i, C_s, cl_msk.sum()), file=prfl)
            # Combine all the masks using a logical OR
            C_masks = np.logical_or.reduce([C_masks, cl_msk])
            N_survived += 1

    return C_masks, N_survived, KDE_vals


def rejctVal(xy, clRjctMethod, KDE_vals, Kest):
    """
    Test how similar this cluster's (x, y) distribution is compared
    to a uniform random distribution.
    """

    if clRjctMethod == 'rkfunc':
        # Test how similar this cluster's (x, y) distribution is compared
        # to a uniform random distribution using Ripley's K.
        # https://stats.stackexchange.com/a/122816/10416

        # Avoid large memory consumption if the data array is too big
        if xy.shape[0] > 5000:
            mode = "none"
        else:
            mode = 'translation'
        rad = np.linspace(.01, .25, 50)

        # Hide RunTimeWarning
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            L_t = Kest.Lfunction(xy, rad, mode=mode)

        # Catch all-nans
        if np.isnan(L_t).all():
            C_s = -np.inf
        else:
            C_s = np.nanmax(abs(L_t - rad))

    elif clRjctMethod == "kdetest":
        from rpy2.robjects import r

        N = xy.shape[0]
        xmin, ymin = xy.min(0)
        xmax, ymax = xy.max(0)

        # Read stored value from table.
        try:
            mean, std = KDE_vals[N]
        except KeyError:
            xrng, yrng = xmax - xmin, ymax - ymin
            r.assign('xrng', xrng)
            r.assign('yrng', yrng)
            r.assign('nstars', N)

            r("""
            maxDistStats <- vector("double", nruns)
            for(i in 1:nruns) {
              dataX <- runif(nstars, 0, xrng)
              dataY <- runif(nstars, 0, yrng)
              kde2dmap <- kde2d(dataX, dataY, n=nKde, lims=c(0, xrng, 0, yrng))
              maxDistStats[i] <- ((max(as.vector(kde2dmap$z))-
              mean(as.vector(kde2dmap$z)))/sd(as.vector(kde2dmap$z)))
            }
            """)
            dist_u = np.array(list(r('maxDistStats')))
            mean, std = np.mean(dist_u), np.std(dist_u)
            KDE_vals[N] = mean, std

        rx = r.matrix(xy.T[0])
        ry = r.matrix(xy.T[1])
        r.assign('minX', xmin)
        r.assign('maxX', xmax)
        r.assign('minY', ymin)
        r.assign('maxY', ymax)
        r.assign('dataX', rx)
        r.assign('dataY', ry)
        r("""
        kde2dmap <- kde2d(dataX, dataY, n=nKde,lims=c(minX, maxX, minY, maxY))
        dist_d <- ((max(as.vector(kde2dmap$z))-mean(as.vector(kde2dmap$z)))/
        sd(as.vector(kde2dmap$z)))
        """)
        dist_d = r('dist_d')[0]

        C_s = (dist_d - mean) / std

    elif clRjctMethod == "kdetestpy":
        N = xy.shape[0]
        xmin, ymin = 0., 0.
        xmax, ymax = xy.max(0)
        xx, yy = np.mgrid[xmin:xmax:50j, ymin:ymax:50j]
        positions = np.vstack([xx.ravel(), yy.ravel()])

        # Read stored value from table.
        try:
            mean, std = KDE_vals[xy.shape[0]]
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
            KDE_vals[N] = mean, std

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

    return C_s, KDE_vals


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
