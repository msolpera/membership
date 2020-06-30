
import numpy as np
from scipy.stats import multivariate_normal


def GUMMtrain(GUMM_perc, clust_xy, probs, prfl, n_epochs=1000, stable_per=.1):
    """
    Fit a model composed of a 2D Gaussian and a 2D uniform distribution in a
    square region with [0., 1.] range.

    Based on the GMM model implementation shown in:
    https://towardsdatascience.com/gaussian-mixture-models-explained-6986aaf5a95
    """
    cl_probs = getProbs(clust_xy, n_epochs, stable_per)

    # Don't overwrite 'probs'
    probs_GUMM = probs.copy()

    # Select the probability cut.
    if GUMM_perc == 'auto':
        from kneebow.rotor import Rotor
        # Create the percentiles (/100.) vs provabilities array.
        percentiles = np.arange(.01, .99, .01)
        perc_probs = np.array([
            percentiles, np.percentile(cl_probs, percentiles * 100.)]).T
        # Find 'knee' where the probabilities start climbing from near 0.
        rotor = Rotor()
        rotor.fit_rotate(perc_probs)
        # Adding 5% to the probability selected here improves the results by
        # making the cut less strict.
        GUMM_prob = perc_probs.T[1][rotor.get_elbow_index()] + 0.05

    else:
        GUMM_prob = np.percentile(cl_probs, GUMM_perc)

    j = 0
    for i, p in enumerate(probs_GUMM):
        # If this was marked as a cluster star
        if p == 1.:
            # And its GUMM probability is below this threshold
            if cl_probs[j] <= GUMM_prob:
                # Mark star as non-member
                probs_GUMM[i] = 0.
            # else:
            #     probs_GUMM[i] = cl_probs[j]
            j += 1

    return probs_GUMM


def getProbs(xy, n_epochs, stable_per):
    """
    """
    cluster = initialize_cluster(xy)

    lkl_old, nstable = -np.inf, 0
    for i in range(n_epochs):

        expectation_step(xy, cluster)
        maximization_step(xy, cluster)

        likelihood = cluster['likelihood']

        # Convergence check 1%
        if abs(likelihood - lkl_old) / likelihood < .1:
            nstable += 1
        if likelihood > lkl_old:
            lkl_old = likelihood
        if nstable == int(stable_per * n_epochs):
            # Converged. Breaking
            break

    # Extract probabilities associated to the 2D Gaussian
    cl_probs = list(cluster['gamma_g'].flatten())

    return cl_probs


def initialize_cluster(data):
    """
    Initialize the 2D Gaussian parameters, and the weights for both
    distributions.
    """
    mu = np.random.uniform(.1, .9, (2,))
    cov = np.eye(2) * np.random.uniform(.1, .9, (2, 2))

    cluster = {'pi_u': .5, 'pi_g': .5, 'mu': mu, 'cov': cov}

    return cluster


def expectation_step(X, cluster):
    """
    """
    # Evaluate Gaussian distribution
    gamma_g = cluster['pi_g'] * multivariate_normal(
        mean=cluster['mu'], cov=cluster['cov']).pdf(X)

    # Evaluate uniform distribution (just a constant)
    gamma_u = cluster['pi_u'] * np.ones(X.shape[0])

    # Normalizing constant
    gammas_sum = gamma_g + gamma_u

    # Probabilities for each element
    cluster['gamma_g'] = gamma_g / gammas_sum
    cluster['gamma_u'] = gamma_u / gammas_sum

    # Save for breaking out
    cluster['likelihood'] = np.sum(np.log(gammas_sum))


def maximization_step(X, cluster):
    """
    """
    gamma_g = cluster['gamma_g']
    N_k = gamma_g.sum(0)

    # Mean
    mu = (gamma_g[:, np.newaxis] * X).sum(0) / N_k
    # Covariance
    cov = np.zeros((X.shape[1], X.shape[1]))
    for j in range(X.shape[0]):
        diff = (X[j] - mu).reshape(-1, 1)
        cov += gamma_g[j] * np.dot(diff, diff.T)
    cov /= N_k

    # Weight for the Gaussian distribution
    N = float(X.shape[0])
    pi_g = N_k / N

    # Weight for the uniform distribution
    pi_u = np.sum(cluster['gamma_u'], axis=0) / N

    # Update parameters
    cluster.update({'pi_u': pi_u, 'pi_g': pi_g, 'mu': mu, 'cov': cov})
