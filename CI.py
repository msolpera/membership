import numpy as np
from scipy.spatial import distance
from scipy import spatial
from astropy.table import Table

def read_data(file_name):
    data = Table.read(file_name, format='ascii')
    data = data['ID', 'x', 'y', 'V', 'BV', 'pmRA', 'pmDE']
    return data['ID'], data['x'], data['y'], data['V'],\
        data['BV'], data['pmRA'], data['pmDE']

def main(rad, cent, file_name):
    # Read data
    ID, coord_x, coord_y, V, BV, pmRA, pmDE = read_data(file_name)

    # Select stars from cluster area (A) and field area (B)
    coords = np.array((coord_x, coord_y)).T
    dist_cent = distance.cdist(cent, coords, 'euclidean')[0]

    dim_a_member = []
    dim_b_member = []
    for i in range(len(coord_x)):
        if dist_cent[i] <= rad:
            dim_a_member.append(pmRA[i])
        else:
            dim_b_member.append(pmRA[i])
    
    # For each star in A, find the 10 nearest neighbors in B and select the one with the greatest distance
    # For each star in A, find the 10 nearest neighbors in A and select de one with the greatest distance
    dist_a_b, dist_a_a = [], []
    rad_ab, rad_aa = [], []
    for i in range(len(dim_a_member)):
        dist_a_b.append(np.abs(dim_a_member[i] - dim_b_member))
        dist_a_a.append(np.abs(dim_a_member[i] - dim_a_member))
    for i in range(len(dim_a_member)):
        rad_ab.append(sorted(dist_a_b[i])[9])
        rad_aa.append(sorted(dist_a_a[i])[10])

    # Calculate the average neighbor density associated with each star of A in B, and idem for A in A
    d_ab, d_aa = [], []
    for i in range(len(dim_a_member)):
        d_ab.append(10/(np.pi*rad_ab[i]**2))
        d_aa.append(10/(np.pi*rad_aa[i]**2))
    d_ab = np.mean(d_ab)
    d_aa = np.mean(d_aa)

    CI = d_ab/d_aa



  
    
    A_parameter, B_paramater = [], []
   

if __name__ == '__main__':
    main(250., [(1000., 1000.)], 'input/synth_clust_out_0.6.dat')
