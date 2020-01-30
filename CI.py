import numpy as np
from scipy.spatial import distance
# from scipy import spatial
from astropy.table import Table
import matplotlib.pyplot as plt


def main(rad, cent, file_name, NN=10):
    # Read data
    data = read_data(file_name)

    # Select stars from cluster area (A) and field area (B)
    coord_x, coord_y = data['x'], data['y']
    coords = np.array((coord_x, coord_y)).T
    dist_cent = distance.cdist(cent, coords, 'euclidean')[0]

    data_dims = ('V', 'BV', 'pmRA', 'pmDE')
    data_arrs = []
    for dim in data_dims:

        # The for loop below can be made faster using a mask, as follows:
        msk = dist_cent <= rad
        dim_a_member = data[dim][msk]
        # The ~ character negates a boolean, i.e., True becomes False
        dim_b_member = data[dim][~msk]

        # dim_a_member = []
        # dim_b_member = []
        # for i in range(len(coord_x)):
        #     if dist_cent[i] <= rad:
        #         dim_a_member.append(pmRA[i])
        #     else:
        #         dim_b_member.append(pmRA[i])

        # The for loop below can be made faster using numpy broadcasting,
        # as follows:
        dist_a_b = np.abs(dim_a_member[:, np.newaxis] - dim_b_member)
        dist_a_a = np.abs(dim_a_member[:, np.newaxis] - dim_a_member)

        # For each star in A, find the NN nearest neighbors in B and select the
        # one with the greatest distance
        # For each star in A, find the NN nearest neighbors in A and select the
        # one with the greatest distance
        # dist_a_b, dist_a_a = [], []
        # for i in range(len(dim_a_member)):
        #     dist_a_b.append(np.abs(dim_a_member[i] - dim_b_member))
        #     dist_a_a.append(np.abs(dim_a_member[i] - dim_a_member))

        # Sort arrays 'in place', so that we don't have to sort them repeatedly
        # inside the for loop below.
        dist_a_b.sort()
        dist_a_a.sort()

        # Calculate the average neighbor density associated with each star of A
        # in B, and idem for A in A
        d_ab, d_aa = [], []
        for i in range(len(dim_a_member)):
            # Distances to the NN neighbor
            rad_ab = dist_a_b[i][NN - 1]
            rad_aa = dist_a_a[i][NN]

            # Store local densities for this star in both regions
            d_ab.append(NN / (np.pi * rad_ab**2))
            d_aa.append(NN / (np.pi * rad_aa**2))

        # Use weighted average, where the weights are the distances to the
        # mean of the B region values (field region)
        weights = abs(np.array(dim_a_member) - np.mean(dim_b_member))
        d_ab_m = np.average(d_ab, weights=weights)
        d_aa_m = np.average(d_aa, weights=weights)

        # Final CI. Divide by the number of stars within the cluster region
        # to normalize.
        CI = (d_ab_m / d_aa_m) / len(d_ab)
        print("{}_CI = {:.2f}".format(dim, CI))

        # Store for plotting
        data_arrs.append([dim_a_member, dim_b_member, CI])

    makePlot(file_name, data_dims, data_arrs)

    # Unused?
    A_parameter, B_paramater = [], []


def read_data(file_name):
    data = Table.read(file_name, format='ascii')
    # data = data['ID', 'x', 'y', 'V', 'BV', 'pmRA', 'pmDE']
    # return data['ID'], data['x'], data['y'], data['V'],\
    #     data['BV'], data['pmRA'], data['pmDE']
    return data


def makePlot(file_name, data_dims, data_arrs):
    """
    """
    fig = plt.figure(figsize=(10, 10))
    for i, (arrA, arrB, CI) in enumerate(data_arrs):

        ax = int("22" + str(i + 1))
        plt.subplot(ax)
        plt.title("CI={:.2f}".format(CI))
        plt.hist(arrA, 25, alpha=.5, color='r', density=True,
                 label="r<rad")
        plt.hist(arrB, 25, alpha=.5, color='b', density=True,
                 label="r>rad")
        plt.xlabel(data_dims[i])
        plt.legend()

    fig.tight_layout()
    CI_txt = file_name[-7:-4]
    plt.savefig("CI_analysis_" + CI_txt + ".png", dpi=150, bbox_inches='tight')


if __name__ == '__main__':
    main(250., [(1024., 1024.)], 'input/synth_clust_out_0.2.dat')
