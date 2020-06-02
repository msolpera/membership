
from astropy.io import ascii
import matplotlib.pyplot as plt
import numpy as np


def main():
    """
    """

    plt.subplot(231)
    resultsProcess(ascii.read("output/up-RESULTS_R_50_1.dat"))
    plt.title("R 1, 50")
    plt.subplot(232)
    resultsProcess(ascii.read("output/up-RESULTS_R_50_2.dat"))
    plt.title("R 2, 50")
    # plt.subplot(233)
    # resultsProcess(ascii.read("output/up-RESULTS_Subl_50_1.dat"))
    # plt.title("Subl 1, 50")
    # plt.subplot(234)
    # resultsProcess(ascii.read("output/up-RESULTS_Subl_50_2.dat"))
    # plt.title("Subl 2, 50")
    plt.subplot(233)
    name = "output/pca_check/up-RESULTS_010_pca2_1.dat"
    print(name)
    resultsProcess(ascii.read(name))
    plt.title("N_PCA=2, 1")
    plt.subplot(234)
    name = "output/pca_check/up-RESULTS_010_pca2_2.dat"
    print(name)
    resultsProcess(ascii.read(name))
    plt.title("N_PCA=2, 2")
    plt.subplot(235)
    name = "output/pca_check/up-RESULTS_010_pca2_3.dat"
    print(name)
    resultsProcess(ascii.read(name))
    plt.title("N_PCA=2, 3")
    # plt.subplot(246)
    # resultsProcess(ascii.read("output/0.1up_50_50.dat"))
    # plt.title("R Sol (N_PCA=4), 50")
    plt.subplot(236)
    resultsProcess(ascii.read("output/up-RESULTS_xycut.dat"))
    plt.title("R xycut, 50")


    # plt.subplot(241)
    # resultsProcess(
    #     ascii.read("output/K_KDE_test/up-0.80_0.49_0.49_0.33_0.29_R_50_1.dat"), ID_col=True)
    # plt.title("R 1, 50")
    # plt.subplot(242)
    # resultsProcess(
    #     ascii.read("output/K_KDE_test/up-0.80_0.49_0.49_0.33_0.29_Subl_50_1.dat"),
    #     ID_col=True)
    # plt.title("R 2, 50")
    # plt.subplot(243)
    # resultsProcess(
    #     ascii.read("output/K_KDE_test/0.80_0.49_0.49_0.33_0.29_R_R.dat"), ID_col=True)
    # plt.title("Python R R, 50")
    # plt.subplot(244)
    # resultsProcess(
    #     ascii.read("output/K_KDE_test/0.80_0.49_0.49_0.33_0.29_R_PR.dat"), ID_col=True)
    # plt.title("Python R PR, 50")
    # plt.subplot(245)
    # resultsProcess(
    #     ascii.read("output/K_KDE_test/0.80_0.49_0.49_0.33_0.29_R_P.dat"), ID_col=True)
    # plt.title("Python R P, 50")
    # plt.subplot(246)
    # resultsProcess(
    #     ascii.read("output/K_KDE_test/0.80_0.49_0.49_0.33_0.29_P_R.dat"), ID_col=True)
    # plt.title("Python P R, 50")
    # plt.subplot(247)
    # resultsProcess(
    #     ascii.read("output/K_KDE_test/0.80_0.49_0.49_0.33_0.29_P_PR.dat"), ID_col=True)
    # plt.title("Python P PR, 50")
    # plt.subplot(248)
    # resultsProcess(
    #     ascii.read("output/K_KDE_test/0.80_0.49_0.49_0.33_0.29_P_P.dat"), ID_col=True)
    # plt.title("Python P P, 50")


    # plt.subplot(231)
    # name = "output/vor_gumm/0.80_0.49_0.49_0.33_0.29_K.dat"
    # print(name)
    # resultsProcess(ascii.read(name), ID_col=True)
    # plt.title("K test, no GUMM")
    # plt.subplot(232)
    # name = "output/vor_gumm/0.80_0.49_0.49_0.33_0.29_KDEpy.dat"
    # print(name)
    # resultsProcess(ascii.read(name), ID_col=True)
    # plt.title("kdetestpy, no GUMM")
    # plt.subplot(233)
    # name = "output/vor_gumm/0.80_0.49_0.49_0.33_0.29_KDER.dat"
    # print(name)
    # resultsProcess(ascii.read(name), ID_col=True)
    # plt.title("kdetest, no GUMM")
    # plt.subplot(234)
    # name = "output/vor_gumm/0.80_0.49_0.49_0.33_0.29_K_G.dat"
    # print(name)
    # resultsProcess(ascii.read(name), ID_col=True)
    # plt.title("K test, GUMM")
    # plt.subplot(235)
    # name = "output/vor_gumm/0.80_0.49_0.49_0.33_0.29_KDEpy_G.dat"
    # print(name)
    # resultsProcess(ascii.read(name), ID_col=True)
    # plt.title("kdetestpy, GUMM")
    # plt.subplot(236)
    # name = "output/vor_gumm/0.80_0.49_0.49_0.33_0.29_KDER_G.dat"
    # print(name)
    # resultsProcess(ascii.read(name), ID_col=True)
    # plt.title("kdetest, GUMM")


    # plt.subplot(231)
    # name = "output/pca_check/up-RESULTS_010_pca2_1.dat"
    # print(name)
    # resultsProcess(ascii.read(name))
    # plt.title("N_PCA=2, 1")
    # plt.subplot(232)
    # name = "output/pca_check/up-RESULTS_010_pca2_2.dat"
    # print(name)
    # resultsProcess(ascii.read(name))
    # plt.title("N_PCA=2, 2")
    # plt.subplot(233)
    # name = "output/pca_check/up-RESULTS_010_pca2_3.dat"
    # print(name)
    # resultsProcess(ascii.read(name))
    # plt.title("N_PCA=2, 3")
    # plt.subplot(234)
    # name = "output/pca_check/0.1up_50_50_ndims4_1.dat"
    # print(name)
    # resultsProcess(ascii.read(name))
    # plt.title("N_PCA=4, 1")
    # plt.subplot(235)
    # name = "output/pca_check/0.1up_50_50_ndims4_2.dat"
    # print(name)
    # resultsProcess(ascii.read(name))
    # plt.title("N_PCA=4, 2")
    # plt.subplot(236)
    # name = "output/pca_check/0.1up_50_50_ndims4_3.dat"
    # print(name)
    # resultsProcess(ascii.read(name))
    # plt.title("N_PCA=4, 3")

    plt.show()


def getMetrics(ID, memb_prob, MP_cut=.9, ID_col=False):
    """
    Obtain the Completeness (C), Purity (P), and Misclassification (M) using
    a cut on MP of 'MP_cut'.
    Also obtain the Logarithmic Scoring Rule (log_MI) which requires no cut.
    """
    # Avoid numerical errors with 0. and 1.
    memb_prob = np.clip(memb_prob, a_min=0.01, a_max=0.99)

    log_MI = 0.
    N_true, N_missed, N_selected, N_interlopers, N_field =\
        0., 0., 0., 0., 0.
    for i in range(len(memb_prob)):

        ID_flag = False
        if ID_col:
            if str(ID[i])[0] == '1':
                ID_flag = True
        else:
            # We are using the 'field' column here, hence the 0
            if str(ID[i]) == '0':
                ID_flag = True

        if ID_flag:
            N_true += 1.
            log_MI += np.log(memb_prob[i])
            if memb_prob[i] >= MP_cut:
                N_selected += 1
            else:
                N_missed += 1
        else:
            N_field += 1
            log_MI += np.log(1. - memb_prob[i])
            if memb_prob[i] >= MP_cut:
                N_interlopers += 1

    # Make it so that 1 is the maximum (best) value, and shorten the range
    # with the logarithm
    L = 1. - np.log(1. - log_MI)

    # Completeness
    C = (N_true - N_missed) / N_true  # = N_selected / N_true
    # Purity
    if N_selected > 0:
        P = (N_selected - N_interlopers) / N_selected
    else:
        P = -np.inf
    # Misclassification
    M = 1. - (N_missed + N_interlopers) / N_true

    print("{} -> C={:.3f}, P={:.3f}, M={:.3f}, L={:.3f}".format(
        MP_cut, C, P, M, L))

    return C, P, M, L, N_true, N_field


def resultsProcess(data, Nruns=50, ID_col=False):

    if ID_col:
        ID_name, prob0, probs_c, prob_f = 'ID', 'prob0', 'prob', 'probs_final'
    else:
        ID_name, prob0, probs_c, prob_f = 'field', 'class', 'class.',\
            'probability'

    metrics = []
    # Store first class
    probs = [data[prob0].data]
    print("{} -> C={:.3f}, P={:.3f}, M={:.3f}, L={:.3f}".format(
        0, *getMetrics(data[ID_name], probs[0], ID_col=ID_col)[:4]))
    for c in range(1, Nruns):
        probs.append(data[probs_c + str(c)].data)
        pr_mean = np.mean(probs, 0)
        C, P, M, L = getMetrics(data[ID_name], pr_mean, ID_col=ID_col)[:4]
        print("{} -> C={:.3f}, P={:.3f}, M={:.3f}, L={:.3f}".format(
            c, C, P, M, L))
        metrics.append([C, P, M, L])

    metrics_final = getMetrics(
        data[ID_name], data[prob_f], ID_col=ID_col)[:4]
    print("Final -> C={:.3f}, P={:.3f}, M={:.3f}, L={:.3f}\n".format(
        *metrics_final))

    metrics = np.array(metrics)
    plt.plot(metrics.T[0], label='C={:.3f}'.format(metrics_final[0]))
    plt.plot(metrics.T[1], label='P={:.3f}'.format(metrics_final[1]))
    plt.plot(metrics.T[2], label='M={:.3f}'.format(metrics_final[2]))
    plt.plot(metrics.T[3], label='L={:.3f}'.format(metrics_final[3]))
    plt.ylim(-2., 1.)
    plt.legend()


if __name__ == '__main__':
    main()
