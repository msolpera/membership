
import numpy as np
from astropy.stats import RipleysKEstimator
import matplotlib.pyplot as plt
import time as t


def RKfunc(Kest, N, rad, mode, N_C_ran=100):
    # https://stats.stackexchange.com/a/122816/10416
    C_vals = []
    for _ in range(N_C_ran):
        xy_u = np.random.uniform(0., 1., (N, 2))
        C_vals.append(Kest(xy_u, (rad,), mode=mode))

    return np.mean(C_vals), np.std(C_vals)


def MCCmetric(TP, FP, TN, FN):
    """
    """
    sqrt = np.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    sqrt = 1. if sqrt == 0. else sqrt
    MCC = ((TP * TN) - (FP * FN)) / sqrt

    return MCC


def successTest():
    """
    Test Ripley's K function to see which configuration has the most
    "separating power" to discriminate between a uniform 2D random
    distribution, and a 2D cluster.
    """

    # If this mode is 'True', it will compare each method in 'K_methods'
    # with the value for a Poisson homogeneous process (which depends on the
    # radius value exclusively).
    poisson_mode = False

    # This is required since under 'poisson_mode' there is no standard
    # deviation and hence the only threshold that means anything is zero.
    if poisson_mode:
        std_list = (0.,)
    else:
        std_list = (1., 2., 3.)

    rad_vals = np.arange(.1, .501, .1)
    N_vals = (10, 25, 50)
    # 'var-width',
    K_methods = ('none', 'ohser', 'ripley', 'translation')

    Kest = RipleysKEstimator(area=1, x_max=1, y_max=1, x_min=0, y_min=0)

    labels = ["{:.2f}".format(_) for _ in rad_vals]
    x = np.arange(len(labels))  # the label locations
    width = 0.2  # the width of the bars
    if poisson_mode:
        widths = (0.,)
    else:
        widths = (-width, 0., width)

    N_run = 2000 if poisson_mode else 100

    for K_method in K_methods:
        print("Method: ", K_method)

        fig, ax123 = plt.subplots(3, 1, figsize=(10, 12))
        start_t = t.time()
        for iN, N in enumerate(N_vals):
            print(" ", N)

            ax123[iN].set_title("N={}".format(N))

            summary = []
            for C_thresh in std_list:
                print("  ", C_thresh)

                C_thresh_lst = []
                for rad in rad_vals:
                    print("   ", rad)

                    # Mean and stddev for a number of RUDs
                    if poisson_mode:
                        mean_k, std_k = Kest.poisson(rad), 1.

                    TP, FP, TN, FN = 0., 0., 0., 0.
                    for _ in range(N_run):

                        if poisson_mode is False:
                            # Process random uniform distributions
                            mean_k, std_k = RKfunc(Kest, N, rad, K_method)

                        # Test a real cluster
                        xy_std = np.random.uniform(0.05, 0.3)
                        xy = np.random.normal((.5, .5), xy_std, (N, 2))
                        # When the observed K value is larger than the
                        # expected K value for a particular distance, the
                        # distribution is more clustered than a random
                        # distribution.
                        C_s = (
                            Kest(xy, (rad,), mode=K_method) - mean_k) / std_k
                        # If C_s is larger than the threshold, the method
                        # correctly identified a real cluster.
                        if C_s >= C_thresh:
                            # Test was successful
                            TP += 1
                        else:
                            # Test failed
                            FN += 1

                        # # Test a RUD
                        # xy_u = np.random.uniform(0., 1., (N, 2))
                        # C_s = (
                        #     Kest(xy_u, (rad,), mode=K_method) - mean_k) / std_k
                        # if C_s > C_thresh:
                        #     FP += 1
                        # else:
                        #     TN += 1

                    # C_thresh_lst.append(MCCmetric(TP, FP, TN, FN))
                    C_thresh_lst.append(TP - FN)

                # summary.append(C_thresh_lst)
                summary.append(np.array(C_thresh_lst) / N_run)

            def autolabel(ax, rects):
                for rect in rects:
                    height = rect.get_height()
                    ax.annotate(
                        '{:.2f}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3), textcoords="offset points", ha='center',
                        va='bottom')

            for iC, Cm in enumerate(summary):
                rect = ax123[iN].bar(
                    x + widths[iC], Cm, width,
                    label=r"$C_s=${}".format(std_list[iC]))
                autolabel(ax123[iN], rect)

            ax123[iN].set_xticks(x)
            ax123[iN].set_xticklabels(labels)
            ax123[iN].set_ylim(.0, 1.1)
            ax123[iN].axhline(1., ls=':', c='k')
            if iN == 0:
                ax123[iN].legend()
            plt.xlabel("rad")
            # ax123[iN].set_ylabel("MCC")
            ax123[iN].set_ylabel("TP-FN")

        elapsed = t.time() - start_t
        plt.suptitle(
            "Method: {}, ({:.0f} secs)".format(K_method, elapsed), y=1.01)
        fig.tight_layout()
        pp = "_poisson" if poisson_mode else ""
        plt.savefig(
            "K_test_{}{}.png".format(K_method, pp), dpi=150,
            bbox_inches='tight')


successTest()
