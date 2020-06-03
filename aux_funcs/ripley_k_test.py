
import numpy as np
from astropy.stats import RipleysKEstimator
import matplotlib.pyplot as plt


def stdRadTest():
    """
    Test Ripley's K function discrimination power using the mean and standard
    deviation of a set of random uniform fields, for several values of the
    'C_s' parameter (distance in standard deviations)
    """

    def RKfunc(Kest, N, rad, N_C_ran=100):
        # https://stats.stackexchange.com/a/122816/10416
        C_vals = []
        for _ in range(N_C_ran):
            xy_u = np.random.uniform(0., 1., (N, 2))
            C_vals.append(Kest(xy_u, (rad,)))

        return np.mean(C_vals), np.std(C_vals)

    # Cluster number of stars (fixed for all runs)
    N = 20
    std_list = (1., 2., 3.)
    rad_lst = (0.5,)
    # Define Ripley's K
    Kest = RipleysKEstimator(area=1, x_max=1, y_max=1, x_min=0, y_min=0)

    results = []
    for C_thresh in std_list:
        print(C_thresh)
        for rad in rad_lst:
            print(rad)
            vals_C_r = []
            for _ in range(100):
                # Standard deviation of the defined cluster
                xy_std = np.random.uniform(0.01, 0.5)
                xy = np.random.normal((.5, .5), xy_std, (N, 2))

                # Process random uniform distributions
                mean_k, std_k = RKfunc(Kest, N, rad)

                # When the observed K value is larger than the expected K
                # value for a particular distance, the distribution is more
                # clustered than a random distribution.
                C_s = (Kest(xy, (rad,)) - mean_k) / std_k

                # We consider 0.25 the limit between cluster and field
                # This is a cluster
                if xy_std <= 0.25:
                    if C_s < C_thresh:
                        # Test failed
                        vals_C_r.append(0.)
                    else:
                        # Test was successful
                        vals_C_r.append(1.)
                # THis is a random field
                else:
                    if C_s < C_thresh:
                        vals_C_r.append(1.)
                    else:
                        vals_C_r.append(0.)

            results.append([
                C_thresh, rad, np.mean(vals_C_r), np.std(vals_C_r)])

    results = np.array(results)
    print(results)


def successTest():
    """
    Test Ripley's K function to see which radius has the most "separating
    power" to discriminate between a uniform 2D random distribution, and a
    2D cluster.
    """

    K_method = 'none'  # 'poisson' 'ohser'  'var-width' 'ripley' 'translation'
    print("Method: ", K_method)
    Kest = RipleysKEstimator(area=1, x_max=1, y_max=1, x_min=0, y_min=0)
    rad_vals = np.linspace(.1, .9, 100)

    for N in (5, 25, 50):
        print(N)
        success = []
        for i, rad in enumerate(rad_vals):
            if not (i + 1) % 50:
                print(" r=", rad)

            if K_method == 'poisson':
                K_un = Kest.poisson(rad)

            rad_res = []
            for _ in range(100):

                # Create a cluster
                xy_std = np.random.uniform(0.05, 0.3)
                xy = np.random.normal((.5, .5), xy_std, (N, 2))

                # Test cluster with a given radius
                K_cl = Kest(xy, (rad,))

                if K_method != 'poisson':
                    # Create random uniform distribution (RUD)
                    xy_u = np.random.uniform(0., 1., (N, 2))
                    # Test RUD with the same radius
                    K_un = Kest(xy_u, (rad,))

                # If it assigned a larger K to the cluster, it means that it
                # was able to identify this as a non RUD; in this case store
                # a value of 1. Otherwise (was not able to correctly identify
                # the cluster) store a value of 0.
                if K_cl > K_un:
                    rad_res.append(1.)
                else:
                    rad_res.append(0.)

            success.append(rad_res)

        # Mean success per radius value.
        success = np.mean(success, 1)
        plt.plot(rad_vals, success, label=r"$N=${}".format(N))

    plt.title("Method: {}".format(K_method))
    plt.legend()
    plt.xlabel("rad")
    plt.ylabel("Success rate")
    plt.show()


# successTest()
stdRadTest()
