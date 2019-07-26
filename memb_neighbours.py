import scipy as sc
from astropy.table import Table
from sklearn import preprocessing
# import pandas as pd
import numpy as np
import seaborn as seaborn
# import matplotlib.pyplot as plt

probando
# Read data
data = Table.read('haf14_match.dat', format='ascii')

# Save data as numpy array
name, coord_x, coord_y, RA, DE, Plx, pmRA, pmDE = np.array(
    [data['id'], data['x'], data['y'], data['RA_ICRS'], data['DE_ICRS'],
     data['Plx'], data['pmRA'], data['pmDE']])

# Normalize input data
RA_norm = preprocessing.normalize(RA[:, np.newaxis], axis=0).ravel()
DE_norm = preprocessing.normalize(DE[:, np.newaxis], axis=0).ravel()
Plx_norm = preprocessing.normalize(Plx[:, np.newaxis], axis=0).ravel()
pmRA_norm = preprocessing.normalize(pmRA[:, np.newaxis], axis=0).ravel()
pmDE_norm = preprocessing.normalize(pmDE[:, np.newaxis], axis=0).ravel()

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
hist_n.set(xlim=(0, 0.02))
hist_n.figure.savefig('hist.png')
