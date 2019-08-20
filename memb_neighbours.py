import scipy as sc
from astropy.table import Table
import pandas as pd
import numpy as np
import seaborn as seaborn
import matplotlib.pyplot as plt
import sys
import time as t

s0 = t.time()
# Read data
data = Table.read('haf14_match.dat', format='ascii')

# Save data as numpy array
name, coord_x, coord_y, RA, DE, pmRA, pmDE = np.array(
    [data['id'], data['x'], data['y'], data['RA_ICRS'], data['DE_ICRS'],
     data['pmRA'], data['pmDE']])
print("1", t.time() - s0)


# Normalize input data


def norm_data(data_array):
    """
    Normalize imput data
    """
    min_array = np.nanmin(data_array)
    max_array = np.nanmax(data_array)
    data_norm = (data_array - min_array)/(max_array - min_array)
    return data_norm


s0 = t.time()
RA_norm = norm_data(RA)
DE_norm = norm_data(DE)
pmRA_norm = norm_data(pmRA)
pmDE_norm = norm_data(pmDE)
print("2", t.time()-s0)

# Assign to each star a vector with its parameters
s0 = t.time()
star = np.array([RA_norm, DE_norm, pmRA_norm, pmDE_norm]).T
print("3", t.time() - s0)

# Calculate the average distance of each star to its nearest n neighbors
s0 = t.time()
star_count = np.zeros(len(name))
n_i = 5
n_f = 25
for n in range(n_i, n_f + 1):
        tree = sc.spatial.cKDTree(star)
        inx = tree.query(star, k=n+1)
        dist = inx[0]
        star_n_prom = []
        for i in range(len(name)):
                promedio = sum(dist[i])/n
                star_n_prom.append(promedio)

        # Stars with distances less than a percentage value are assigned the  
        # number 1 (Cluster member)
        p_i = 1
        p_f = 10
        for p in range(p_i, p_f + 1):
                d_max = np.percentile(star_n_prom, p)
                for j in range(len(name)):
                        if star_n_prom[j] < d_max:
                                star_count[j] = star_count[j] + 1
print("4", t.time() - s0)

# A percentage of membership is assigned to each star
star_count_max = np.max(star_count)
star_m_perc = star_count/star_count_max

# The coordinates of the stars with percentage of membership greater than perc 
# are plotted
memb_RA = []
memb_DE = []
perc = 0.5
for i in range(len(name)):
        if star_m_perc[i] > perc:
                memb_DE.append(DE[i])
                memb_RA.append(RA[i])

plt.scatter(memb_RA, memb_DE)
plt.show()
