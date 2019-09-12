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
    memb_prob = np.zeros(n_stars)
    n_i, n_f = 5, 25
    for n in range(n_i, n_f + 1):
        star_n_prom = dist_neighbours(n, star, n_stars)

        # Calculate the median of the average distances
        median = np.median(star_n_prom)

        # Select stars with average distance less than median
        select_star = []
        for i in range(n_stars):
            if star_n_prom[i] < median:
                select_star.append(star[i])
        
        # Calculate the 4th dimensional centroid of these stars
        centroid = sum(select_star)/len(select_star)
        
        # Calculate the probability of membership as the reciprocal of the
        # distance to the centroid
        for i in range(n_stars):
            cent_dist = np.linalg.norm(star[i] - centroid)
            memb_prob[i] = memb_prob[i] + 1/cent_dist
        
    # Calculate the average of probability of membership and normalize
    mp_prom = np.array(memb_prob/(n_f-n_i))
    mp_norm = norm_data(mp_prom)
    
    # Generate plot
    plot_memb(mp_norm, DE, RA, n_stars)


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


def plot_memb(mp_norm, DE, RA, n_stars):
    memb_RA, memb_DE, memb_color = [], [], []
    for i in range(n_stars):
        if mp_norm[i] > 0.7:
            memb_DE.append(DE[i])
            memb_RA.append(RA[i])
            memb_color.append(mp_norm[i])

    plt.title("Number of stars: N={}".format(len(memb_RA)))
    plt.scatter(memb_RA, memb_DE, c=memb_color)
    plt.colorbar()
    plt.savefig('membership_out.png', dpi=150, bbox_inches='tight')


if __name__ == '__main__':
    main()
