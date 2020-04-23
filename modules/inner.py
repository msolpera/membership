
import numpy as np
from sklearn.cluster import KMeans, SpectralClustering, MiniBatchKMeans,\
    AgglomerativeClustering, DBSCAN, OPTICS
from sklearn.mixture import GaussianMixture
from scipy.cluster.vq import kmeans2
from .voronoiVols import voronoi_volumes
from . import compFunc


def main(
    clust_xy, clust_data, clust_method, otlrFlag, C_thresh, unif_method,
        RK_rad, clust_params, vol_cummul):
    """
    """

    # Obtain all the clusters in the input data using kMeans
    clusts_msk = clustAlgor(clust_data, clust_method, clust_params)

    # Reject "field" clusters and return their stored masks
    C_masks = rjctField(
        unif_method, RK_rad, vol_cummul, C_thresh, clust_xy, clusts_msk)

    return C_masks, len(clusts_msk)


def clustAlgor(clust_data, clust_method, clust_params):
    """
    Find 'n_clusters' in the 'data' array using kMeans.
    """

    # Number of clusters: min is 2, max is N_cl_max
    n_clusters = max(2, min(
        int(clust_data.shape[0] / float(clust_params['N_membs'])),
        int(clust_params['N_cl_max'])))

    if clust_method in ('kmeans', 'minibatchKMeans'):

        kinit_method = 'random'

        if kinit_method == 'voronoi_1':
            # Select the centers after clustering the most dense cells. Not
            # good performance apparently
            data_vols = voronoi_volumes(clust_data[:, :3], False)
            msk = data_vols < np.median(data_vols)
            if n_clusters > msk.sum():
                centroid, n_clusters = clust_data[msk], msk.sum()
            else:
                centroid, _ = kmeans2(clust_data[msk], n_clusters)
            cents, n_init = centroid, 1

        elif kinit_method == 'voronoi_2':
            # This selects the points with the smallest Voronoi volumes as the
            # centers. Bad purity
            data_vols = voronoi_volumes(clust_data[:, :3], False)
            idx = np.argpartition(data_vols, n_clusters)[:n_clusters]
            cents, n_init = clust_data[idx], 1

        elif kinit_method in ('k-means++', 'random'):
            cents, n_init = kinit_method, 10

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
            if probs and len(probs) == clust_data.shape[0]:
                print("iterating probs")
                data_vols = voronoi_volumes(clust_data[:, :3], False)
                vals, weights = [], []
                for i, p in enumerate(probs):
                    if p > .5:
                        vals.append(clust_data[i])
                        weights.append(1. / data_vols[i])

                k_cent = np.average(vals, axis=0, weights=weights)

                # Assign rest of centers randomly
                centers = np.random.uniform(
                    clust_data.min(0), clust_data.max(0),
                    [n_clusters - 1, clust_data.shape[1]])

                # Assign rest of centers using Voronoi (worse P)
                # idx = np.argpartition(
                #    data_vols, n_clusters - 1)[:n_clusters - 1]
                # centers = data[idx]

                cents, n_init = np.concatenate([centers, [k_cent]]), 1
            else:
                cents, n_init = 'k-means++', n_init

        if clust_method == 'kmeans':
            # KMeans
            model = KMeans(
                n_clusters=n_clusters, init=cents, n_init=n_init, verbose=0)
        else:
            model = MiniBatchKMeans(
                n_clusters=n_clusters, init=cents, n_init=n_init, verbose=0)

    elif clust_method == 'spectral':
        model = SpectralClustering(n_clusters=n_clusters)

    elif clust_method == 'agglomerative':
        # Faster than kMeans, deterministic (?) --> better C but lower P
        model = AgglomerativeClustering(n_clusters=n_clusters)

    elif clust_method == 'GMM':
        model = GaussianMixture(n_components=n_clusters)

    elif clust_method == 'DBSCAN':
        # # Finding eps manually:
        # # https://towardsdatascience.com/
        # # machine-learning-clustering-dbscan-determine-the-optimal-value-for-
        # # epsilon-eps-python-example-3100091cfbc
        # Another method described in:
        # Amin Karami and Ronnie Johansson. Choosing dbscan parameters
        # automatically using differential evolution. International Journal
        # of Computer Applications, 91(7), 2014
        #
        # from sklearn.neighbors import NearestNeighbors
        # neigh = NearestNeighbors(n_neighbors=2)
        # nbrs = neigh.fit(data)
        # distances, indices = nbrs.kneighbors(data)
        # distances = np.sort(distances, axis=0)
        # distances = distances[:,1]
        # import matplotlib.pyplot as plt
        # plt.plot(distances)
        # # Finding eps with 'kneed':
        # # https://github.com/arvkevi/kneed
        model = DBSCAN(eps=.2)

    elif clust_method == 'OPTICS':
        model = OPTICS()

    model.fit(clust_data)
    if clust_method == 'GMM':
        labels = model.fit_predict(clust_data)
    else:
        labels = model.labels_
    # print(model.inertia_)

    # Scipy's implementation of kmeans. It is much slower
    # centroid, labels = kmeans2(data, n_clusters, iter=max_iter)

    # Separate the labels that point to each cluster found
    clusts_msk = [(labels == _) for _ in range(labels.max() + 1)]

    print(" N stars={}, N clusters={}".format(
        clust_data.shape[0], len(clusts_msk)))

    return clusts_msk


def rjctField(unif_method, RK_rad, vol_cummul, C_thresh, clust_xy, clusts_msk):
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
            # the areas of a a uniform random distribution.
            C_s = compFunc.main(
                unif_method, RK_rad, clust_xy[cl_msk], vol_cummul, vol_d)

        # Smaller C_s values point to samples that come from the
        # same distribution (i.e., this is a "field cluster")
        if C_s < C_thresh:
            # Store mask that points to stars that should be *removed*
            # (hence we "flip" it with the '~' operator)
            C_masks.append(~cl_msk)

    return C_masks
