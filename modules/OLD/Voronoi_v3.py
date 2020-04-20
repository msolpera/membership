
import warnings
import numpy as np

from scipy.stats import gaussian_kde

from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KernelDensity, KNeighborsClassifier
from sklearn.mixture import GaussianMixture
from sklearn.naive_bayes import GaussianNB, BernoulliNB
from sklearn import svm
from sklearn.linear_model import LogisticRegression
from sklearn.cluster import KMeans, SpectralClustering, AgglomerativeClustering
from sklearn import metrics

from voronoiVols import voronoi_volumes
from dataRead import read_data
from update_progress import updt

import matplotlib.pyplot as plt


def main(file_name):
    """
    """

    # Read data, normalized without outliers.
    ID, stars = read_data(file_name)

    # Calculate Voronoi volumes
    vol = voronoi_volumes(stars)

    # vol_median = np.median(vol)
    # msk = vol < vol_median
    # stars = stars[msk]
    # vol = vol[msk]

    # # Select the optimal percentile that separates the member and
    # # field stars.
    # # First grid search
    # obj_func_val, opt_vol, optm_perc = gridSearch(
    #     stars, vol, 0., 0., 0., ext=50., N_grid=50)
    # # Second grid search
    # obj_func_val, opt_vol, optm_perc = gridSearch(
    #     stars, vol, obj_func_val, opt_vol, optm_perc, ext=5.)
    # # Third (final) grid search
    # _, opt_vol, optm_perc = gridSearch(
    #     stars, vol, obj_func_val, opt_vol, optm_perc, ext=1.)
    # print(opt_vol, optm_perc)

    if file_name.startswith('input/0.1'):
        opt_vol, optm_perc = 0.0021182989083568266, 34.10646258503401
    if file_name.startswith('input/0.9'):
        opt_vol, optm_perc = 6.154918184405654e-06, 0.7801764455782313
        # opt_vol, optm_perc = 6.112082676368036e-06, 1.4715136054421771
# 
    for method in ('kNN', 'NB', 'BB', 'SVM', 'logit'):
        print("Method: {}".format(method))
        # Initial field priors.
        field_priors = np.ones(stars.shape[0])
        # KDE evaluation
        p_m, p_f = evalKDE(stars, vol, optm_perc, field_priors, method, run_1rst=True)
        # Bayesian probabilities
        _, field_priors = BayesProbs(vol, field_priors, p_m, p_f, opt_vol)

        # Run several times re-sampling the data given its uncertainties?
        for _ in np.ones(10) * .75:  #np.linspace(0.75, 0.5, 10):
            # Re-sample input data (re-use Voronoi volumes)
            # stars = reSample(stars)

            p_m, p_f = evalKDE(stars, vol, optm_perc, field_priors, method, p_cut=_)
            probability, field_priors = BayesProbs(
                vol, field_priors, p_m, p_f, opt_vol)

            from process_synth_clust import member_index
            print("P={:.2f}, MI={:.5f}, TI_90={:.5f}".format(
                _, *member_index(ID, probability)[:2]))

        # Plot stars with probability greater than 0.9
        x_memb, y_memb, pmra_memb, pmde_memb, V_memb, BV_memb, prob = \
            [], [], [], [], [], [], []
        for l, st in enumerate(stars):
            if probability[l] > 0.5:
                x_memb.append(st[0])
                y_memb.append(st[1])
                pmra_memb.append(st[2])
                pmde_memb.append(st[3])
                # V_memb.append(V[l])
                # BV_memb.append(BV[l])
                prob.append(probability[l])
        print(len(x_memb))
        plt.subplot(221)
        plt.hist(probability) # , density=True)
        plt.subplot(222)
        plt.scatter(x_memb, y_memb, s=20., c=prob, lw=.5, edgecolor='k')
        plt.colorbar(aspect=90, pad=0.01)
        plt.subplot(223)
        plt.scatter(pmra_memb, pmde_memb, s=20., c=prob, lw=.5, edgecolor='k')
        plt.colorbar(aspect=90, pad=0.01)
        # plt.subplot(224)
        # plt.scatter(BV_memb, V_memb, s=20., c=prob, lw=.5, edgecolor='k')
        # plt.gca().invert_yaxis()
        # plt.colorbar(aspect=90, pad=0.01)
        plt.show()

    return ID, probability


def gridSearch(stars, vol, obj_func_val, opt_vol, optm_perc, ext, N_grid=25):
    """
    """
    pmin, pmax = max(0.05, optm_perc - ext), min(50., optm_perc + ext)

    for i, perc in enumerate(np.linspace(pmin, pmax, N_grid)):
        obj_func_val_n, opt_vol_n = objFunc(stars, vol, perc)
        if obj_func_val_n > obj_func_val:
            obj_func_val, opt_vol, optm_perc = obj_func_val_n, opt_vol_n, perc
        updt(N_grid, i + 1)

    print('Perc, selected, objFval: {:.2f}, {:.2f}, {:.2f}'.format(
        perc, optm_perc, obj_func_val))

    return obj_func_val, opt_vol, optm_perc


def objFunc(stars, vol, perc):
    """
    Objective function to select the optimal percentile.
    """
    vol_p = np.percentile(vol, perc)
    msk_memb = vol <= vol_p

    # Stars below the volume percentile
    X = stars[msk_memb]

    obj_val = 0.
    if X.shape[0] < 10:
        # # The 'davies_bouldin_score' minimizes to 0.
        # obj_val = 10.
        pass
    else:
        with warnings.catch_warnings():
            # warnings.filterwarnings('ignore')
            warnings.simplefilter("ignore")

            # # KMeans
            # model = KMeans(
            #     n_clusters=2, init='k-means++', random_state=1, n_init=100)

            # SpectralClustering
            model = SpectralClustering(
                n_clusters=2, affinity='nearest_neighbors', n_init=100,
                random_state=1, assign_labels='discretize')

            # # Hierarchical clustering
            # model = AgglomerativeClustering(n_clusters=2, linkage='ward')

            # # # Gaussian MM
            # from sklearn import mixture
            # model = mixture.GaussianMixture(
            #     n_components=2, covariance_type='full')
            # model.fit(X)
            # probs = model.predict_proba(stars)

        model.fit(X)
        labels = model.labels_

        try:
            # More indexes listed here:
            # https://stats.stackexchange.com/a/79080/10416

            # Evaluate clustering
            obj_val = metrics.silhouette_score(X, labels, metric='euclidean')
            # obj_val = metrics.calinski_harabasz_score(X, labels)
            # obj_val = -metrics.davies_bouldin_score(X, labels)
        except ValueError:
            # 'silhouette_score' (possibly other too?) throw this exception if
            # the number of labels is 1.
            pass

    return obj_val, vol_p


def evalKDE(
    stars, vol, optm_perc, field_priors, method, run_1rst=False,
        p_cut=.5):
    """
    p_cut :
    """

    # Separate stars given the optimal percentile
    vol_p = np.percentile(vol, optm_perc)
    msk_memb = vol <= vol_p
    membs_stars = stars[msk_memb]
    field_stars = stars[~msk_memb]

    # # This fails for very contaminated clusters that can not be clearly
    # # separated using Voronoi because the mean is dominated by field
    # # stars.
    #
    # memb_cent, memb_std = membs_stars.mean(0), membs_stars.std()
    # from scipy.spatial import distance
    # dd = distance.cdist([memb_cent], membs_stars)
    # msk_dd = dd[0] < 2. * memb_std
    # membs_stars = membs_stars[msk_dd]

    # Skip first run
    if not run_1rst:
        priors_membs = 1. - field_priors[msk_memb]
        qrt_1_membs = priors_membs < p_cut

        # If there are any stars classified as members by the Voronoi mask,
        # with probabilities (estimated as: 1. - field_prob) less than the
        # minimum allowed value (i.e., 'p_cut'), then make them field stars.
        memb_prior = 1. - field_priors
        if qrt_1_membs.sum() > 1:
            for i, flag in enumerate(msk_memb):
                if flag:
                    # Purported cluster star
                    if memb_prior[i] < p_cut:
                        # Flip flag to field star
                        msk_memb[i] = False

        priors_field = field_priors[~msk_memb]
        qrt_1_field = priors_field < p_cut

        if qrt_1_field.sum() > 1:
            for i, flag in enumerate(~msk_memb):
                if flag:
                    # Purported field star
                    if field_priors[i] < p_cut:
                        # Flip flag to cluster member
                        msk_memb[i] = True

    # Separate stars using the new mask
    membs_stars = stars[msk_memb]
    field_stars = stars[~msk_memb]

    if membs_stars.shape[0] < 2 or field_stars.shape[0] < 2:
        return np.ones(stars.shape[0]) * .5, np.ones(stars.shape[0]) * .5

    # if method == 'KDE':
    #     Nst_max=1000
    #     # To improve the performance, cap the number of stars using a random
    #     # selection of 'Nf_max' elements.
    #     if field_stars.shape[0] > Nst_max:
    #         idxs = np.arange(field_stars.shape[0])
    #         np.random.shuffle(idxs)
    #         field_stars = field_stars[idxs[:Nst_max]]
    #     # Evaluate all stars in both KDEs
    #     try:
    #         bw_f = sklearnBW(field_stars)
    #         kd_field = gaussian_kde(field_stars.T, bw_method=bw_f)
    #         bw_m = sklearnBW(membs_stars)
    #         kd_memb = gaussian_kde(membs_stars.T, bw_method=bw_m)
    #         # TODO Bottleneck is here
    #         p_memb = kd_memb.evaluate(stars.T)
    #         p_field = kd_field.evaluate(stars.T)
    #     except (np.linalg.LinAlgError, ValueError):
    #         return np.ones(stars.shape[0]), np.ones(stars.shape[0])
    #     # Normalize to make the probabilities behave better, particularly for
    #     # large contaminations.
    #     # WHY DOES THIS WORK?
    #     p_memb /= p_memb.max()
    #     p_field /= p_field.max()

    if method == 'kNN':
        # KNeighborsClassifier
        model = KNeighborsClassifier()  # weights="distance"
    elif method == 'NB':
        # Naive Bayes
        model = GaussianNB()
    elif method == 'BB':
        # # Bernoulli Bayes
        model = BernoulliNB()
    elif method == 'SVM':
        # SVM. This methods requires probability=True
        model = svm.SVC(probability=True) #class_weight={1: 10})
    elif method == 'logit':
        # Logit regression
        model = LogisticRegression()

    # vol_membs = vol[msk_memb]
    # vol_field = vol[~msk_memb]
    # def volW(x):
    #     w_vols = np.concatenate([vol_membs, vol_field])
    #     w_vols = np.array([w_vols, w_vols, w_vols, w_vols, w_vols]).T
    #     return w_vols

    # Create a 'training set' with the stars separated by the above process.
    X = np.concatenate([membs_stars, field_stars])
    # Create hard labels for both sets.
    labels = [0 for _ in range(membs_stars.shape[0])] +\
        [1 for _ in range(field_stars.shape[0])]
    # Train the model using the training set
    model.fit(X, labels)

    # print("Score: {:.3f}".format(model.score(X, labels)))

    # This requires checking which cluster corresponds to the actual cluster
    # in a post-processing step.
    # # fit a Gaussian Mixture Model with two components
    # model = GaussianMixture(n_components=2, covariance_type='full')
    # # Train the model using the training set
    # model.fit(X)

    try:
        # Predict probabilities.
        probs = model.predict_proba(stars)
        p_memb, p_field = probs.T
    except:
        import pdb; pdb.set_trace()  # breakpoint 0f500a05 //


    # plt.subplot(121)
    # plt.scatter(membs_stars.T[0], membs_stars.T[1], c='g', zorder=5)
    # plt.scatter(field_stars.T[0], field_stars.T[1], c='orange')
    # plt.subplot(122)
    # plt.scatter(membs_stars.T[2], membs_stars.T[3], c='g', zorder=5)
    # plt.scatter(field_stars.T[2], field_stars.T[3], c='orange')

    return p_memb, p_field


def sklearnBW(data, rtol=1.e-3, bwmin=.01, bwmax=.5, bwN=10, cv=20):
    """
    Bandwidth selection.

    http://jakevdp.github.io/blog/2013/12/01/kernel-density-estimation/
    """
    # 20-fold cross-validation
    grid = GridSearchCV(
        KernelDensity(), {'bandwidth': np.linspace(bwmin, bwmax, bwN),
                          'rtol': (rtol,)}, cv=cv)
    grid.fit(data)
    bandwidth = grid.best_params_['bandwidth']
    # Needed for compatibility with scipy's gaussian_kde
    bw = bandwidth / data.std(ddof=1)
    return bw


def BayesProbs(vol, field_priors, p_m, p_f, opt_vol):
    """
    Calculate the Bayesian probability of membership
    """

    # with warnings.catch_warnings():
    #     warnings.filterwarnings('ignore')

    #     # Gaussian volume prior
    #     sigma = opt_vol
    #     prior_prob_cl = (1. / sigma) * np.exp(-.5 * (vol / sigma)**2)

    #     # TODO WHY DOES THIS WORK???
    #     prior_prob_cl = prior_prob_cl / prior_prob_cl.max()

    #     prior = field_priors / prior_prob_cl

    #     # # Flat prior
    #     # prior = 1.

    #     probability = 1. / (1. + prior * (p_f / p_m))
    #     field_priors = 1. - probability

    probability = p_m
    field_priors = p_f

    return probability, field_priors


if __name__ == '__main__':
    main()
