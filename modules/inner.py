
import numpy as np
import sklearn.cluster as skclust
import sklearn.mixture as skmixture
from astropy.stats import RipleysKEstimator


def main(
    clust_xy, clust_data, clust_method, RK_rad, C_thresh, clust_params,
        cl_method_pars):
    """
    Perform the inner loop: cluster --> reject
    """
    print("  Performing clustering on array of shape ({}, {})".format(
        *clust_data.shape))

    # Obtain all the clusters in the input data using kMeans
    clusts_msk = clustAlgor(
        clust_data, clust_method, clust_params, cl_method_pars)
    print("  Identified {} clusters".format(len(clusts_msk)))

    # Reject "field" clusters and return their stored masks
    C_masks = rjctField(clust_xy, clusts_msk, RK_rad, C_thresh)

    return C_masks, len(clusts_msk)


def clustAlgor(clust_data, clust_method, clust_params, cl_method_pars):
    """
    Find 'n_clusters' in the 'data' array using kMeans.
    """

    if clust_method == 'KMeans':
        model = skclust.KMeans()

    elif clust_method == 'MiniBatchKMeans':
        model = skclust.MiniBatchKMeans()

    elif clust_method == 'AffinityPropagation':
        model = skclust.AffinityPropagation()

    elif clust_method == 'SpectralClustering':
        model = skclust.SpectralClustering()

    elif clust_method == 'AgglomerativeClustering':
        model = skclust.AgglomerativeClustering()

    elif clust_method == 'GaussianMixture':
        model = skmixture.GaussianMixture()

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
        model = skclust.DBSCAN()

    elif clust_method == 'OPTICS':
        model = skclust.OPTICS()

    elif clust_method == 'MeanShift':
        model = skclust.MeanShift()

    elif clust_method == 'Birch':
        model = skclust.Birch()

    elif clust_method == 'Voronoi':
        # print("Obtaining Voronoi volumes")
        from .voronoiVols import voronoi_volumes
        from scipy.spatial.distance import cdist
        from kneed import KneeLocator
        vol_v = voronoi_volumes(clust_data, False)
        dens = 1. / vol_v
        dist = cdist(clust_data, clust_data)
        delta = np.zeros(dens.size)
        for i, st_dens in enumerate(dens):
            msk = dens > st_dens
            if msk.sum() == 0:
                idx_max = i
            else:
                delta[i] = dist[i][msk].min()
        delta[idx_max] = delta.max()
        mult = dens * delta
        # Indexes that sort in descending order
        idx_s = np.argsort(-mult)

    # Set parameters for the method (if any)
    if cl_method_pars:
        model.set_params(**cl_method_pars)

    # Only these methods require the number of clusters to be set
    if clust_method in (
            'KMeans', 'MiniBatchKMeans', 'SpectralClustering',
            'AgglomerativeClustering', 'GaussianMixture'):

        # Number of clusters: min is 2, max is N_cl_max
        n_clusters = max(2, min(
            int(clust_data.shape[0] / clust_params['N_membs']),
            clust_params['N_cl_max']))

        if clust_method == 'GaussianMixture':
            model.n_components = n_clusters
        else:
            model.n_clusters = n_clusters

    elif clust_method == 'Voronoi':
        # Locate the 'knee' in the curve, to estimate the number of clusters
        kneedle = KneeLocator(
            range(len(mult)), mult[idx_s], S=10, curve='convex',
            direction='decreasing')
        n_clusters = max(2, min(int(kneedle.knee), clust_params['N_cl_max']))

    if clust_method == 'Voronoi':
        # Obtain labels
        labels = np.argmin(dist[idx_s[:n_clusters], :], 0)
    else:
        # Fit the model
        model.fit(clust_data)

        # Extract the labels
        if clust_method == 'GaussianMixture':
            labels = model.fit_predict(clust_data)
        else:
            labels = model.labels_

    # Separate the labels that point to each cluster found
    clusts_msk = [(labels == _) for _ in range(labels.min(), labels.max() + 1)]

    return clusts_msk


def rjctField(clust_xy, clusts_msk, RK_rad, C_thresh, N_C_ran=100):
    """
    Iterate over all the clusters defined by the clustering algorithm, storing
    those masks that point to clusters that survived the test with unniform
    random 2D fields in (x, y.
    """

    # Test how similar this cluster's (x, y) distribution is compared
    # to a uniform random distribution using Ripley's K.
    # https://stats.stackexchange.com/a/122816/10416
    # Define the (x, y) area with sides [0, 1]
    Kest = RipleysKEstimator(area=1, x_max=1, y_max=1, x_min=0, y_min=0)

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
            C_vals = []
            for _ in range(N_C_ran):
                xy_u = np.random.uniform(0., 1., clust_xy[cl_msk].shape)
                C_vals.append(Kest(xy_u, (RK_rad,)))

            # When the observed K value is larger than the expected K value for
            # a particular distance, the distribution is more clustered than a
            # random distribution at that distance
            mean, std = np.mean(C_vals), np.std(C_vals)
            C_s = (Kest(clust_xy[cl_msk], (RK_rad,)) - mean) / std

        # Smaller C_s values point to samples that come from a uniform
        # distribution in (x, y), i.e., a "cluster" made up of field
        # stars. Hence, we keep as "true" clusters those with C_s values
        # larger than C_thresh.
        if C_s >= C_thresh:
            # Store mask that points to stars that should be *kept*
            print("   Cluster {} survived (C_s={:.2f}), N={}".format(
                i, C_s[0], cl_msk.sum()))
            C_masks.append(cl_msk)

    return C_masks
