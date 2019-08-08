import scipy as sc
from astropy.table import Table
import pandas as pd
import numpy as np
import seaborn as seaborn
import matplotlib.pyplot as plt
import sys

# Read data
data = Table.read('haf14_match.dat', format='ascii')

# Save data as numpy array
name, coord_x, coord_y, RA, DE, Plx, pmRA, pmDE = np.array(
    [data['id'], data['x'], data['y'], data['RA_ICRS'], data['DE_ICRS'],
     data['Plx'], data['pmRA'], data['pmDE']])


# msk = (pmRA < 0) & (pmRA > -4) & (pmDE < 5) & (pmDE > 0) 
# df = pd.DataFrame(np.array([RA, DE]).T, columns=["RA_ICRS", "DE_ICRS"])
# df2 =  pd.DataFrame(np.array([pmRA[msk],pmDE[msk]]).T, columns=["pmRA", "pmDE"])
# plt.scatter(RA, DE)
# seaborn.jointplot(x="RA_ICRS", y="DE_ICRS", data = df, kind="kde")
# seaborn.jointplot(x="pmRA", y="pmDE", data = df2, kind="kde",)
# plt.show()
# plt.scatter(pmRA, pmDE)
# plt.show()
# sys.exit()


# Normalize input data

def norm_data(data_array):
    """
    Normalize imput data
    """
    min_array = np.nanmin(data_array)
    max_array = np.nanmax(data_array)
    data_norm = (data_array - min_array)/(max_array - min_array)
    return data_norm

RA_norm = norm_data(RA)
DE_norm = norm_data(DE)
Plx_norm = norm_data(Plx)
pmRA_norm = norm_data(pmRA)
pmDE_norm = norm_data(pmDE)


# Assign to each star a vector with its parameters

star = np.array([RA_norm, DE_norm, Plx_norm, pmRA_norm, pmDE_norm]).T

# star = []
# for i in range(len(name)):
# star.append((RA_norm[i], DE_norm[i],
# Plx_norm[i], pmRA_norm[i], pmDE_norm[i]))


# Calculate distance between stars

dist = sc.spatial.distance.cdist(star, star, 'euclidean')


# Sort distances from least to greatest

dist_star_sort = np.sort(dist)
# dist_star = []
# for i in range(len(name)):
#    dist_star.append(dist[i])
#    dist_star[i].sort()


# Select the n smallest distances

n = 10
dist_star_n_min = []
for i in range(len(name)):
    min = dist_star_sort[i]
    min_n = min[1:n + 1]
    dist_star_n_min.append(min_n)


# Average of the n shortest distances

star_n_prom = []
for i in range(len(name)):
    prom = sum(dist_star_n_min[i]) / n
    star_n_prom.append(prom)

# Assign a number from 0 to 1 to each star (1 = min dist)

star_n_prom_norm = np.array(norm_data(star_n_prom))
star_n_prom_factor = star_n_prom_norm*(-1.)+1.

memb_RA = []
memb_DE = []
for i in range(len(name)):
    if star_n_prom_factor[i] >= 0.9:
        memb_RA.append(RA[i])
        memb_DE.append(DE[i])

print(len(memb_RA))

# plt.subplot(121)
# plt.hist(star_n_prom, bins=50)
# plt.subplot(122)
# plt.hist(star_n_prom_factor, bins=50)
plot = seaborn.JointGrid(memb_RA, memb_DE, space=0, size=5, ratio=30)
plot.plot_joint(plt.scatter, color="darkcyan")
plt.show()
