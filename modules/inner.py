
import numpy as np
from sklearn.cluster import KMeans, SpectralClustering, MiniBatchKMeans,\
    AgglomerativeClustering
from scipy.cluster.vq import kmeans2
from .voronoiVols import voronoi_volumes
from . import compFunc


def main(
    clust_xy, cl_data, clust_method, kM_N_membs, kM_N_cl_max, kM_n_init,
        kM_max_iter, N_C_ran, method, vol_cummul, C_thresh, probs):
    """
    """

    # Obtain all the clusters in the input data using kMeans
    clusts_msk = clustAlgor(
        cl_data, clust_method, kM_N_membs, kM_N_cl_max, kM_n_init, kM_max_iter,
        probs)

    # Reject "field" clusters and return their stored masks
    C_masks = rjctField(
        N_C_ran, method, vol_cummul, C_thresh, clust_xy, clusts_msk)

    return C_masks, len(clusts_msk)


def clustAlgor(
        data, clust_method, N_membs, N_cl_max, n_init, max_iter, probs):
    """
    Find 'n_clusters' in the 'data' array using kMeans.
    """

    # Number of clusters: min is 2, max is N_cl_max
    n_clusters = max(2, min(int(data.shape[0] / N_membs), N_cl_max))

    if clust_method in ('kmeans', 'minibatchKMeans'):

        kinit_method = 'k-means++'

        if kinit_method == 'voronoi_1':
            # Select the centers after clustering the most dense cells. Not
            # good performance apparently
            data_vols = voronoi_volumes(data[:, :3], False)
            msk = data_vols < np.median(data_vols)
            if n_clusters > msk.sum():
                centroid, n_clusters = data[msk], msk.sum()
            else:
                centroid, _ = kmeans2(data[msk], n_clusters, iter=max_iter)
            cents, n_init = centroid, 1

        elif kinit_method == 'voronoi_2':
            # This selects the points with the smallest Voronoi volumes as the
            # centers. Bad purity
            data_vols = voronoi_volumes(data[:, :3], False)
            idx = np.argpartition(data_vols, n_clusters)[:n_clusters]
            cents, n_init = data[idx], 1

        elif kinit_method == 'k-means++':
            cents, n_init = 'k-means++', n_init

        elif kinit_method == 'iterative':
            # Run after the first outer loop, each first iteration of the inner
            # loop. The idea is to use the stars identified as members in the
            # previous iteration to estimate the center of the cluster, using
            # random centers for the rest of the 'n_clusters'.
            # It is not clear if it works or not. E.g., for the synth clust
            # with CI=0.8 the k-means++ method gives:
            # C=0.595, P=0.000, log_MI=-538
            # and this one:
            # C=0.608, P=-0.067, log_MI=-480
            if probs and len(probs) == data.shape[0]:
                print("iterating probs")
                data_vols = voronoi_volumes(data[:, :3], False)
                vals, weights = [], []
                for i, p in enumerate(probs):
                    if p > .5:
                        vals.append(data[i])
                        weights.append(1. / data_vols[i])

                k_cent = np.average(vals, axis=0, weights=weights)

                # Assign rest of centers randomly
                centers = np.random.uniform(
                    data.min(0), data.max(0), [n_clusters - 1, data.shape[1]])

                # Assign rest of centers using Voronoi (worse P)
                # idx = np.argpartition(data_vols, n_clusters - 1)[:n_clusters - 1]
                # centers = data[idx]

                cents, n_init = np.concatenate([centers, [k_cent]]), 1
            else:
                cents, n_init = 'k-means++', n_init

        if clust_method == 'kmeans':
            # KMeans
            model = KMeans(
                n_clusters=n_clusters, init=cents, n_init=n_init,
                max_iter=max_iter, verbose=0)
        else:
            model = MiniBatchKMeans(
                n_clusters=n_clusters, init=cents, n_init=n_init,
                max_iter=max_iter, verbose=0)

    elif clust_method == 'spectral':
        model = SpectralClustering(n_clusters=n_clusters)

    elif clust_method == 'agglomerative':
        # Faster than kMeans, deterministic (?) --> better C but lower P
        model = AgglomerativeClustering(n_clusters=n_clusters)

    model.fit(data)
    labels = model.labels_
    # print(model.inertia_)

    # Scipy's implementation of kmeans. It is much slower
    # centroid, labels = kmeans2(data, n_clusters, iter=max_iter)

    # Separate the labels that point to each cluster found
    clusts_msk = [(labels == _) for _ in range(labels.max() + 1)]

    print(" N stars={}, N clusters={}".format(
        data.shape[0], len(clusts_msk)))

    return clusts_msk


def rjctField(N_C_ran, method, vol_cummul, C_thresh, clust_xy, clusts_msk):
    """
    Iterate over all the clusters defined by the kMeans algorithm, storing
    those masks that point to clusters rejected as not real clusters.
    """

    C_masks = []
    # For each cluster found by the kMeans, check if it is composed
    # of field stars or actual cluster members, using their (x, y) Voronoi
    # area distribution, and the selected test.
    for i, cl_msk in enumerate(clusts_msk):

        if cl_msk.sum() < 5:
            # HARDCODED
            # If the cluster contains less than 5 stars, skip the test and
            # mark its stars as field stars
            C_s = -np.inf
        else:
            # Obtain the Voronoi areas in (x, y) for the stars in this cluster
            vol_d = voronoi_volumes(clust_xy[cl_msk])

            # Test how similar this cluster's area distribution is compared to
            # the areas of a a uniform random distribution. Repeat for
            # 'N_C_ran' random uniform distributions for better statistics,
            # keeping the minimum value at each loop.
            C_s = compFunc.main(
                method, N_C_ran, clust_xy[cl_msk], vol_cummul, vol_d)

        # Smaller C_s values point to samples that come from the
        # same distribution (i.e., this is a "field cluster")
        if C_s < C_thresh:
            # Store mask that points to stars that should be *removed*
            # (hence we "flip" it with the '~' operator)
            C_masks.append(~cl_msk)

    return C_masks
