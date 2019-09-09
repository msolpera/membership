from scipy import spatial
from astropy.table import Table
# import pandas as pd
import numpy as np
# import seaborn as seaborn
import matplotlib.pyplot as plt
# import sys
import time as t


def main(ID, RA, DE, pmRA, pmDE):

    # Number of stars
    n_stars = len(ID)

    # Normalize input data
    RA_norm = norm_data(RA)
    DE_norm = norm_data(DE)
    pmRA_norm = norm_data(pmRA)
    pmDE_norm = norm_data(pmDE)

    # Assign to each star a vector with its parameters
    star = np.array([RA_norm, DE_norm, pmRA_norm, pmDE_norm]).T

    # Calculate the average distance of each star to its nearest n neighbors
    star_count = np.zeros(n_stars)
    n_i, n_f = 5, 25
    for n in range(n_i, n_f + 1):
        star_n_prom = dist_neighbours(n, star, n_stars)
        # Stars with distances less than a percentage value are assigned the
        # number 1 (Cluster member)
        p_i, p_f = 1, 10
        for p in range(p_i, p_f + 1):
            d_max = np.percentile(star_n_prom, p)
            for j in range(n_stars):
                if star_n_prom[j] < d_max:
                    star_count[j] = star_count[j] + 1

    # A percentage of membership is assigned to each star
    star_count_max = np.max(star_count)
    star_m_perc = star_count / star_count_max

    # The coordinates of the stars with percentage of membership greater than
    # perc are plotted
    perc = 0.75
    plot_memb(perc, DE, RA, star_m_perc, n_stars)


def norm_data(data_array):
    """
    Normalize input data
    """
    min_array = np.nanmin(data_array)
    max_array = np.nanmax(data_array)
    data_norm = (data_array - min_array) / (max_array - min_array)
    return data_norm


def dist_neighbours(n, star, n_stars):
    tree = spatial.cKDTree(star)
    inx = tree.query(star, k=n)  # <-- This is the correct number
    dist = inx[0]
    star_n_prom = []
    for i in range(n_stars):
        promedio = sum(dist[i]) / n  # <-- np.mean(dist[i]) could be used
        star_n_prom.append(promedio)
    return(star_n_prom)


def plot_memb(perc, DE, RA, star_m_perc, n_stars):
    memb_RA, memb_DE, memb_color = [], [], []
    for i in range(n_stars):
        if star_m_perc[i] > perc:
            memb_DE.append(DE[i])
            memb_RA.append(RA[i])
            memb_color.append(star_m_perc[i])

    plt.title("Number of stars: N={}".format(len(memb_RA)))
    plt.scatter(memb_RA, memb_DE, c=memb_color)
    plt.colorbar()
    plt.savefig('membership_out.png', dpi=150, bbox_inches='tight')



if __name__ == '__main__':
    main()
