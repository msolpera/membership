import scipy as sc
from astropy.table import Table
from sklearn import preprocessing
# import pandas as pd
import numpy as np
import seaborn as seaborn
import matplotlib.pyplot as plt

# Read data
data = Table.read('haf14_match.dat', format='ascii')

# Save data as numpy array
name, coord_x, coord_y, RA, DE, Plx, pmRA, pmDE = np.array(
    [data['id'], data['x'], data['y'], data['RA_ICRS'], data['DE_ICRS'],
     data['Plx'], data['pmRA'], data['pmDE']])

# Normalize input data


def norm_data(data_array):
    """
    Normalize imput data
    """
    min_array = np.nanmin(data_array)
    data_min = data_array - min_array
    max_array = np.nanmax(data_min)
    data_norm = data_min/max_array
    return data_norm

RA_norm = norm_data(RA)
DE_norm = norm_data(DE)
Plx_norm = norm_data(Plx)
pmRA_norm = norm_data(pmRA)
pmDE_norm = norm_data(pmDE)


# Assign to each star a vector with its parameters
l = 1.
star = []
for i in range(len(name)):
    star.append((l * RA_norm[i], l * DE_norm[i],
                 Plx_norm[i], pmRA_norm[i], pmRA_norm[i]))

# Calculate distance between stars
dist = sc.spatial.distance.cdist(star, star, 'euclidean')
# print(dist)

# Sort distances from least to greatest
dist_star = []
for i in range(len(name)):
    dist_star.append(dist[i, :])
    dist_star[i].sort()
# print(dist_star[0])

# Select the n smallest distances

n = 10
dist_star_n_min = []
for i in range(len(name)):
    min = dist_star[i]
    min_n = min[1:n + 1]
    dist_star_n_min.append(min_n)
# print(dist_star_n_min[0])

# Average of the n shortest distances

star_n_prom = []
for i in range(len(name)):
    prom = sum(dist_star_n_min[i]) / n
    star_n_prom.append(prom)

# print(star_n_prom)

hist_n = seaborn.distplot(star_n_prom, bins=100, kde=False, color='g')
plt.show()
# hist_n.figure.savefig('hist.png')
