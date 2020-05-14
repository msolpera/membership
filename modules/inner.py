
import numpy as np
import sklearn.cluster as skclust
import sklearn.mixture as skmixture


def main(
    clust_xy, clust_data, clust_method, RK_rad, RK_vals, Kest, C_thresh,
        clust_params, cl_method_pars, prfl):
    """
    Perform the inner loop: cluster --> reject
    """
    print("  Performing clustering on array of shape ({}, {})".format(
        *clust_data.shape), file=prfl)

    # Obtain all the clusters in the input data using kMeans
    clusts_msk = clustAlgor(
        clust_data, clust_method, clust_params, cl_method_pars)
    print("  Identified {} clusters".format(len(clusts_msk)), file=prfl)

    # Reject "field" clusters and return their stored masks
    C_masks = rjctField(
        clust_xy, clusts_msk, RK_rad, RK_vals, Kest, C_thresh, prfl)

    return C_masks, len(clusts_msk)


def clustAlgor(clust_data, clust_method, clust_params, cl_method_pars):
    """
    Find clusters in the 'clust_data' array using the selected algorithm.
    """

    # Number of clusters: min is 2, max is N_cl_max
    n_clusters = max(2, min(
        int(clust_data.shape[0] / clust_params['N_membs']),
        clust_params['N_cl_max']))

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

    elif clust_method == 'BayesianGaussianMixture':
        model = skmixture.BayesianGaussianMixture()

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
        from sklearn.neighbors import NearestNeighbors
        neigh = NearestNeighbors(n_neighbors=cl_method_pars['min_samples'])
        nbrs = neigh.fit(clust_data)
        distances, indices = nbrs.kneighbors(clust_data)
        distances = np.sort(distances, axis=0)
        distances = distances[:, 1]
        # import matplotlib.pyplot as plt
        # plt.plot(distances)
        # Finding eps with 'kneed'
        from kneed import KneeLocator
        S = cl_method_pars['knee_s']
        kneedle = KneeLocator(
            range(len(distances)), distances, S=S, curve='convex',
            direction='increasing')
        eps = kneedle.knee_y
        if eps is None:
            # Select maximum value as eps
            eps = distances[-1]
        model = skclust.DBSCAN(eps=eps)

    elif clust_method == 'OPTICS':
        model = skclust.OPTICS()

    elif clust_method == 'MeanShift':
        model = skclust.MeanShift()

    elif clust_method == 'Birch':
        model = skclust.Birch()

    elif clust_method == 'Voronoi':
        # Adapted from: 'Clustering by fast search and find of density peaks',
        # Rodriguez and Laio (2014)
        from .voronoiVols import voronoi_volumes
        from scipy.spatial.distance import cdist
        from kneed import KneeLocator
        # Obtain Voronoi volumes
        vol_v = voronoi_volumes(clust_data, False)
        # Convert to densities
        dens = 1. / vol_v
        # Find the distances to all stars, for all stars
        dist = cdist(clust_data, clust_data)
        # For each star, find the distance to the closest star that has a
        # larger density. For the star with largest density, assign the
        # distance to the most distant star.
        delta = np.zeros(dens.size)
        for i, st_dens in enumerate(dens):
            msk = dens > st_dens
            # Store the index of the star with the largest density.
            if msk.sum() == 0:
                idx_max = i
            else:
                delta[i] = dist[i][msk].min()
        # For this star, assign the largest distance.
        delta[idx_max] = delta.max()
        mult = dens * delta
        # Indexes that sort in descending order
        idx_s = np.argsort(-mult)

    elif clust_method == 'HDBSCAN':
        import hdbscan
        model = hdbscan.HDBSCAN()

    elif clust_method == 'KMeansSwap':
        model = skclust.KMeans()
        model.n_clusters = n_clusters
        model.fit(clust_data)

        inertia_old = model.inertia_
        centers = model.cluster_centers_
        for _ in range(cl_method_pars['n_runs']):
            centers2 = np.array(centers)

            idx_1 = np.random.choice(n_clusters)
            idx_2 = np.random.choice(clust_data.shape[0])
            centers2[idx_1] = clust_data[idx_2]

            model = skclust.KMeans(
                init=centers2, n_clusters=n_clusters, n_init=1, max_iter=2)
            model.fit(clust_data)
            if model.inertia_ < inertia_old:
                centers = model.cluster_centers_
                inertia_old = model.inertia_

        # Reset this parameter
        model.max_iter = 300

    elif clust_method == 'pyclKmeans':
        from pyclustering.cluster.kmeans import kmeans
        from pyclustering.cluster.center_initializer import kmeans_plusplus_initializer

        initial_centers = kmeans_plusplus_initializer(
            clust_data, n_clusters).initialize()
        model = kmeans(clust_data, initial_centers)
        # final_centers = model.get_centers()

    elif clust_method == 'pyclGA':
        from pyclustering.cluster.ga import genetic_algorithm
        # Create instance of observer that will collect all information:
        # observer_instance = ga_observer(True, True, True)
        model = genetic_algorithm(
            clust_data, count_clusters=n_clusters, chromosome_count=100,
            population_count=20, coeff_mutation_count=.5)

    #
    # Set parameters for the method (if any)
    if cl_method_pars and clust_method not in ('Voronoi', 'KMeansSwap'):
        if clust_method == 'DBSCAN':
            cl_method_pars2 = {
                k: cl_method_pars[k] for k in cl_method_pars if k != 'knee_s'}
            model.set_params(**cl_method_pars2)
        else:
            model.set_params(**cl_method_pars)

    # Only these methods require the number of clusters to be set
    if clust_method in (
            'KMeans', 'MiniBatchKMeans', 'SpectralClustering',
            'AgglomerativeClustering', 'GaussianMixture',
            'BayesianGaussianMixture'):

        if clust_method in ('GaussianMixture', 'BayesianGaussianMixture'):
            model.n_components = n_clusters
        else:
            model.n_clusters = n_clusters

    elif clust_method == 'Voronoi':
        # Locate the 'knee' in the curve, to estimate the number of clusters
        S = cl_method_pars['knee_s']
        kneedle = KneeLocator(
            range(len(mult)), mult[idx_s], S=S, curve='convex',
            direction='decreasing')
        n_clusters = max(2, min(int(kneedle.knee), clust_params['N_cl_max']))

    # Fit the model
    if clust_method[:4] == 'pycl':
        model.process()
    elif clust_method != 'Voronoi':
        # Fit the model
        model.fit(clust_data)

    # Extract the labels
    if clust_method in ('GaussianMixture', 'BayesianGaussianMixture'):
        labels = model.predict(clust_data)
        # probs_gmm = model.predict_proba(clust_data).max(1)
    elif clust_method == 'Voronoi':
        labels = np.argmin(dist[idx_s[:n_clusters], :], 0)
    elif clust_method[:4] == 'pycl':
        if clust_method[4:] == 'Kmeans':
            labels = model.predict(clust_data)
        elif clust_method[4:] == 'GA':
            labels = np.zeros(clust_data.shape[0])
            for i, clust in enumerate(model.get_clusters()):
                labels[clust] = i
            labels = labels.astype(int)
    else:
        labels = model.labels_

    # Separate the labels that point to each cluster found
    clusts_msk = [(labels == _) for _ in range(labels.min(), labels.max() + 1)]

    # import matplotlib.pyplot as plt
    # for cl in clusts_msk:
    #     plt.scatter(*clust_data[cl].T, alpha=.5)

    return clusts_msk


def rjctField(
    clust_xy, clusts_msk, RK_rad, RK_vals, Kest, C_thresh, prfl,
        minStars=5):
    """
    Iterate over all the clusters defined by the clustering algorithm, storing
    those masks that point to clusters that survived the test with uniform
    random 2D fields in (x, y).

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
            C_s = rejctVal(clust_xy[cl_msk], RK_rad, RK_vals, Kest)

        # Smaller C_s values point to samples that come from a uniform
        # distribution in (x, y), i.e., a "cluster" made up of field
        # stars. Hence, we keep as "true" clusters those with C_s values
        # larger than C_thresh.
        if C_s >= C_thresh:
            # Store mask that points to stars that should be *kept*
            print("   Cluster {} survived (C_s={:.2f}), N={}".format(
                i, C_s, cl_msk.sum()), file=prfl)
            C_masks.append(cl_msk)

    return C_masks


def rejctVal(xy, RK_rad, RK_vals, Kest, N_C_ran=100, method='3'):
    """
    Test how similar this cluster's (x, y) distribution is compared
    to a uniform random distribution using Ripley's K.

    https://stats.stackexchange.com/a/122816/10416
    """

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

    # elif method == '3':
    # This method compares the original xy distribution with lots of
    # uniform random samples
    try:
        mean, std = RK_vals[xy.shape[0]]
    except KeyError:
        # If the number of stars in this cluster is larger than the largest
        # one stored, used the values for this largest value.
        mean, std = list(RK_vals.values())[-1]

    # for _ in range(N_C_ran):
    #     xy_u = np.random.uniform(0., 1., xy.shape)
    #     C_vals.append(Kest(xy_u, (RK_rad,)))
    # # When the observed K value is larger than the expected K value for
    # # a particular distance, the distribution is more clustered than a
    # # random distribution at that distance
    # mean, std = np.mean(C_vals), np.std(C_vals)

    C_s = (Kest(xy, (RK_rad,))[0] - mean) / std

    return C_s
