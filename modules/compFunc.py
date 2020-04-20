
import numpy as np
from astropy.stats import RipleysKEstimator
import warnings
from .voronoiVols import voronoi_CDF2vols, voronoi_vols2CDF


def main(method, N_C_ran, xy, vol_cummul, vol_d):
    """
    Use the selected test to evaluate how similar the Voronoi areas for the
    analyzed cluster are compared to a random uniform distribution.
    """
    if method == 'MV1':
        # "Mean volume" method
        C_vals = []
        for _ in range(N_C_ran):
            # THIS DOES NOT WORK!
            # xy_u = np.random.uniform(xy.min(0), xy.max(0), xy.shape)
            xy_u = np.random.uniform(0., 1., xy.shape)
            vol_u = voronoi_CDF2vols(xy_u, vol_cummul)
            C_vals.append(vol_u.mean())

        mean, std = np.mean(C_vals), np.std(C_vals)
        C_s = (mean - vol_d.mean()) / std

    elif method == 'MV2':
        C_vals = []
        for _ in range(N_C_ran):
            xy_u = np.random.uniform(0., 1., xy.shape)
            vol_u = voronoi_CDF2vols(xy_u, vol_cummul)
            C_vals.append((vol_u.mean() - vol_u.min()) / vol_u.std())

        C_vol_d = (vol_d.mean() - vol_d.min()) / vol_d.std()
        C_s = (C_vol_d - np.mean(C_vals)) / np.std(C_vals)

    elif method == 'AD1':
        # Avoid printing lots of useless warning messages
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            C_vals = []
            for _ in range(N_C_ran):
                # Sample areas for random uniform (x, y) data with matching
                # number of elements.
                # IMPORTANT: This assumes that the (x, y) data is
                # normalized to (0, 1) in both dimensions!
                xy_u = np.random.uniform(0., 1., xy.shape)
                vol_u1 = voronoi_CDF2vols(xy_u, vol_cummul)
                vol_u2 = voronoi_CDF2vols(xy_u, vol_cummul)
                C_u = anderson_ksamp_new([vol_u1, vol_u2])
                C_d1 = anderson_ksamp_new([vol_d, vol_u1])
                C_d2 = anderson_ksamp_new([vol_d, vol_u2])
                C_vals.append([C_d1, C_d2, C_u])

        C_d1, C_d2, C_u = np.array(C_vals).T
        AD_mean, AD_std = np.concatenate([C_d1, C_d2]).mean(),\
            np.concatenate([C_d1, C_d2]).std()
        AD_unif = C_u.mean()
        C_s = (AD_mean - AD_unif) / AD_std

    elif method == 'AD2':
        # Avoid printing lots of useless warning messages
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            C_vals = []
            for _ in range(N_C_ran):
                # Sample areas for random uniform (x, y) data with matching
                # number of elements.
                xy_u = np.random.uniform(0., 1., xy.shape)
                vol_u1 = voronoi_CDF2vols(xy_u, vol_cummul)
                vol_u2 = voronoi_CDF2vols(xy_u, vol_cummul)
                C_u = anderson_ksamp_new([vol_u1, vol_u2])
                C_d1 = anderson_ksamp_new([vol_d, vol_u1])
                C_d2 = anderson_ksamp_new([vol_d, vol_u2])
                C_vals.append([C_d1, C_d2, C_u])

            C_d1, C_d2, C_u = np.array(C_vals).T
            C_s = anderson_ksamp_new([np.concatenate([C_d1, C_d2]), C_u])

    elif method == 'KS':
        # Bootstrap resample of the volumes
        vol_d = np.random.choice(vol_d, vol_d.size)
        C_s = kstest(vol_d, voronoi_vols2CDF, args=(vol_cummul,))

    elif method == 'RK':
        # https://stats.stackexchange.com/a/122816/10416
        Kest = RipleysKEstimator(area=1, x_max=1, y_max=1, x_min=0, y_min=0)
        rad = .5
        C_vals = []
        for _ in range(N_C_ran):
            xy_u = np.random.uniform(0., 1., xy.shape)
            C_vals.append(Kest(xy_u, (rad,)))

        # When the observed K value is larger than the expected K value for
        # a particular distance, the distribution is more clustered than a
        # random distribution at that distance
        mean, std = np.mean(C_vals), np.std(C_vals)
        C_s = (Kest(xy, (rad,)) - mean) / std

    # if C_s > 3.:
    #     # from scipy.spatial.distance import cdist
    #     # from .voronoiVols import voronoi_CDF2vols
    #     from scipy.spatial import Voronoi, voronoi_plot_2d
    #     import matplotlib.pyplot as plt

    #     vor = Voronoi(xy)

    #     plt.subplot(221)
    #     plt.scatter(xy.T[0], xy.T[1], c=vol_d)
    #     plt.colorbar()
    #     ax = plt.subplot(222)
    #     voronoi_plot_2d(vor, ax=ax)
    #     plt.subplot(223)
    #     plt.hist(C_vals, 25)
    #     plt.axvline(x=C_vol_d, c='g', zorder=5)
    #     plt.axvline(x=mean, c='r', zorder=5)
    #     plt.axvline(x=mean - std, ls='--', c='r', zorder=5)

    return C_s


def anderson_ksamp_new(samples):
    """
    Source: https://github.com/scipy/scipy/blob/v1.4.1/scipy/stats/
            morestats.py#L1914-L2070

    A2akN: equation 7 of Scholz and Stephens.

    samples : sequence of 1-D array_like
        Array of sample arrays.
    Z : array_like
        Sorted array of all observations.
    Zstar : array_like
        Sorted array of unique observations.
    k : int
        Number of samples.
    n : array_like
        Number of observations in each sample.
    N : int
        Total number of observations.
    Returns
    -------
    A2aKN : float
        The A2aKN statistics of Scholz and Stephens 1987.

    """

    k = 2  # len(samples)
    # samples = list(map(np.asarray, samples))

    Z = np.sort(np.hstack(samples))
    N = Z.size
    Zstar = np.unique(Z)

    n = np.array([sample.size for sample in samples])

    # midrank = True
    # A2kN = _anderson_ksamp_midrank(samples, Z, Zstar, k, n, N)

    A2kN = 0.
    Z_ssorted_left = Z.searchsorted(Zstar, 'left')
    if N == Zstar.size:
        lj = 1.
    else:
        lj = Z.searchsorted(Zstar, 'right') - Z_ssorted_left
    Bj = Z_ssorted_left + lj / 2.
    for i in np.arange(0, k):
        s = np.sort(samples[i])
        s_ssorted_right = s.searchsorted(Zstar, side='right')
        Mij = s_ssorted_right.astype(float)
        fij = s_ssorted_right - s.searchsorted(Zstar, 'left')
        Mij -= fij / 2.
        inner = lj / float(N) * (N*Mij - Bj*n[i])**2 / (Bj*(N - Bj) - N*lj/4.)
        A2kN += inner.sum() / n[i]
    A2kN *= (N - 1.) / N

    H = (1. / n).sum()
    hs_cs = (1. / np.arange(N - 1, 1, -1)).cumsum()
    h = hs_cs[-1] + 1
    g = (hs_cs / np.arange(2, N)).sum()

    a = (4*g - 6) * (k - 1) + (10 - 6*g)*H
    b = (2*g - 4)*k**2 + 8*h*k + (2*g - 14*h - 4)*H - 8*h + 4*g - 6
    c = (6*h + 2*g - 2)*k**2 + (4*h - 4*g + 6)*k + (2*h - 6)*H + 4*h
    d = (2*h + 6)*k**2 - 4*h*k
    sigmasq = (a*N**3 + b*N**2 + c*N + d) / ((N - 1.) * (N - 2.) * (N - 3.))
    m = k - 1
    A2 = (A2kN - m) / np.sqrt(sigmasq)

    return A2


def kstest(rvs, cdf, args=()):
    """
    """

    vals = np.sort(rvs)
    N = len(vals)
    cdfvals = cdf(vals, *args)

    Dplus = (np.arange(1.0, N + 1) / N - cdfvals).max()
    Dmin = (cdfvals - np.arange(0.0, N) / N).max()

    D = np.max([Dplus, Dmin])
    return D
