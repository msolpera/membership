
import numpy as np
from scipy.spatial import distance
from astropy.table import Table, vstack
from astropy.io import ascii
import matplotlib.pyplot as plt


def main(CI=0.8):
    """

    xy_range : float
      The x,y range of the total frame in pixels.
    cl_area: float
      The area of the cluster region. Set to be a circle of r=250 px.
    CI : float < 1.
      The contamination index. Set to be 0.6 to match Haffner 14.
    """

    # Hardcoded parameters: cluster's center and radius, and frame's lengths
    # in x & y.
    cl_cent, cl_rad, xy_range = [(1024., 1024.)], 250., 2048.

    # Read input synthetic cluster data
    data = ascii.read('synth_clust_input.dat')

    # Extract members and field stars data
    membs_ID, membs_x, membs_y, membs_V, membs_BV, field_V, field_BV =\
        membSeparate(data)
    N_membs = len(membs_ID)

    # cluster area given the radius defined
    cl_area = np.pi * cl_rad**2
    # Total area of the square frame, given the x,y range defined
    tot_area = xy_range**2

    # Estimate number of members given the CI
    N_field = estimateNfield(N_membs, CI, tot_area, cl_area)

    # Generate field stars IDs
    field_ID = generateFielfIDs(N_field)

    # Generate position data for the field stars.
    x_fl, y_fl = np.random.uniform(0., xy_range, (2, N_field))

    # Generate V, BV data for the field stars, using the photometric data
    # from the members as template.
    field_V_N, field_BV_N = generateFieldPhot(N_field, field_V, field_BV)

    # Generate proper motions data
    pms_membs, pms_field, dmean = generatePMs(N_membs, N_field)

    # Write final combined data to file
    tabl_cl = Table(
        [membs_ID, membs_x, membs_y, membs_V, membs_BV, pms_membs[0],
         pms_membs[1]], names=('ID', 'x', 'y', 'V', 'BV', 'pmRA', 'pmDE'))
    tabl_fl = Table(
        [field_ID, x_fl, y_fl, field_V_N, field_BV_N, pms_field[0],
         pms_field[1]], names=('ID', 'x', 'y', 'V', 'BV', 'pmRA', 'pmDE'))
    
    # Generate plot
    makePlot(
        cl_cent, cl_rad, CI, membs_x, membs_y, membs_V, membs_BV, x_fl, y_fl,
        field_BV_N, field_V_N, pms_membs, pms_field, dmean)
    
    # Generate CI for another dimensions
    table_data = vstack([tabl_cl, tabl_fl])
    CI_array, data_dims, data_arrs = CI_2(cl_rad, cl_cent, table_data, NN=10)

    # Gererate plor for another dimensions
    makePlot_2(CI, data_dims, data_arrs)

    # Write final data
    ascii.write(
        vstack([tabl_cl, tabl_fl]), str(CI) + '_' + str(CI_array[0]) + '_' + str(CI_array[1]) 
        + '_' + str(CI_array[2]) + '_' + str(CI_array[3]) + '.dat',
        overwrite=True)

    print("Finished")


def membSeparate(data):
    """
    Separate actual members and field stars by their IDs, assuming real
    members have IDs that start with a '1'.
    """
    membs_ID, membs_x, membs_y, membs_V, membs_BV = [[] for _ in range(5)]
    field_V, field_BV = [[] for _ in range(2)]
    for st in data:
        if str(st['id'])[0] == '1':
            # This is a real member
            membs_ID.append(str(st['id']))
            membs_x.append(st['x'])
            membs_y.append(st['y'])
            membs_V.append(st['V'])
            membs_BV.append(st['BV'])
        else:
            field_V.append(st['V'])
            field_BV.append(st['BV'])

    return membs_ID, membs_x, membs_y, membs_V, membs_BV, field_V, field_BV


def estimateNfield(N_membs, CI, tot_area, cl_area):
    """
    Estimate the total number of field stars that should be generated so
    that the CI is respected.
    """

    # Number of field stars in the cluster area
    N_field_in_clreg = N_membs / ((1. / CI) - 1.)         #CONSULTAR

    # Field stars density
    field_dens = N_field_in_clreg / cl_area

    # Total number of field stars in the entire frame
    N_field = int(field_dens * tot_area)

    return N_field


def generateFielfIDs(N_field):
    """
    Assign IDs for field stars, all starting with '2'.
    """
    field_ID = []
    for i in range(N_field):
        field_ID.append('2' + str(i))

    return field_ID


def generateFieldPhot(N_field, field_V, field_BV):
    """
    Generate synthetic V vs BV data for the field stars, by randomly
    displacing the field stars in the input file.
    """
    field_V_N, field_BV_N = [], []
    for i in range(N_field):
        field_V_N.append(np.random.choice(field_V) + np.random.normal(0., .1))
        field_BV_N.append(
            np.random.choice(field_BV) + np.random.normal(0., .1))

    return field_V_N, field_BV_N


def generatePMs(
    N_membs, N_field, fm1=-5., fm2=5., fstd1=1., fstd2=5., cm1=.05,
        cm2=.1, nstd=3.):
    """
    Given a 'data' array read from a cluster file, generate proper motions
    from a bi-variate Gaussian with reasonable mean and standard deviation
    values. Also generates reasonable uncertainties, based on the values
    taken from the Gaia DR2 survey.

    Parameters
    ----------
    N_membs : int
       Number of members stars.
    N_field : int
       Number of field stars.
    fm1 : float
       Minimum value for the `mean` range used for field stars.
    fm2 : float
       Maximum value for the `mean` range used for field stars.
    fstd1 : float
       Minimum value for the `stddev` range used for field stars.
    fstd2 : float
       Maximum value for the `stddev` range used for field stars.
    cm1 : float
       Minimum value for the scale used for the `stddev` of cluster stars.
    cm2 : float
       Maximum value for the scale used for the `stddev` of cluster stars.
    nstd : float
       Number of field stars standard deviations around the mean, that define
       the range where the mean for the cluster stars can be generated.
    """

    # Generate PMs for field stars.
    mean_f = np.random.uniform(fm1, fm2, (2))
    # This defines the standard deviation for the field stars PMs.
    c1 = np.random.uniform(fstd1, fstd2)
    # Defined as a covariance matrix, since this is a bi-variate Gaussian.
    cov = np.eye(2) * c1
    # Generate PMs for all the field stars.
    pms_field = np.random.multivariate_normal(mean_f, cov, N_field)

    # Generate PMs for member stars.
    c2 = np.random.uniform(cm1, cm2)
    # Use the same covariance matrix used for the field stars, but scaled with
    # a factor 'c2'.
    cov = cov * c2

    # For the cluster members, use a `mean` taken randomly from the limits of
    # the PM values already defined for the field stars.

    # Limits in the field stars PM_RA
    ra_std = pms_field.T[0].std()
    range_RA = [mean_f[0] - nstd * ra_std, mean_f[0] + nstd * ra_std]
    # Random mean in the RA axis, within the defined range
    m1 = np.random.uniform(*range_RA)

    # Limits in the field stars PM_DEC
    de_std = pms_field.T[0].std()
    range_DE = [mean_f[1] - nstd * de_std, mean_f[1] + nstd * de_std]
    # Random mean in the DEC axis, within the defined range
    m2 = np.random.uniform(*range_DE)

    # Mean for cluster stars
    mean_c = np.array([m1, m2])
    # Generate PMs for all the member stars.
    pms_membs = np.random.multivariate_normal(mean_c, cov, N_membs)

    # Distance between means, normalized by the cluster's region stddev.
    dmean = np.linalg.norm(mean_f - mean_c) / cov[0][0]

    return pms_membs.T, pms_field.T, dmean


def makePlot(
    cl_cent, cl_rad, CI, membs_x, membs_y, membs_V, membs_BV, x_fl, y_fl,
        field_BV_N, field_V_N, pms_membs, pms_field, dmean):
    """
    """
    dist_cent = distance.cdist(cl_cent, np.array((x_fl, y_fl)).T)[0]
    msk = dist_cent <= cl_rad

    fig = plt.figure(figsize=(10, 15))
    plt.suptitle("CI={:.2f}, d_mean={:.2f}".format(CI, dmean), y=1.02)

    plt.subplot(321)
    x_fl_i, y_fl_i = x_fl[msk], y_fl[msk]
    x_fl_o, y_fl_o = x_fl[~msk], y_fl[~msk]
    plt.scatter(
        x_fl_o, y_fl_o, c='b', ec='k', lw=.5, alpha=.7, s=10,
        label="Field")
    plt.scatter(
        x_fl_i, y_fl_i, c='orange', ec='k', lw=.5, alpha=.7, s=25,
        label="Field, r<rad")
    plt.scatter(membs_x, membs_y, c='r', s=40, ec='k', label="Cluster")
    plt.legend(fontsize="small")

    plt.subplot(322)
    plt.hexbin(
        x_fl.tolist() + membs_x, y_fl.tolist() + membs_y, gridsize=20,
        cmap='Greens')

    plt.subplot(323)
    pmra_i, pmde_i = pms_field[0][msk], pms_field[1][msk]
    pmra_o, pmde_o = pms_field[0][~msk], pms_field[1][~msk]
    plt.scatter(pmra_o, pmde_o, c='b', s=10, ec='k', lw=.5, alpha=.7,)
    plt.scatter(pmra_i, pmde_i, c='orange', s=25, ec='k', lw=.5, alpha=.7,)
    plt.scatter(pms_membs[0], pms_membs[1], c='r', s=30, ec='k')

    plt.subplot(324)
    plt.hexbin(
        pms_field[0].tolist() + pms_membs[0].tolist(),
        pms_field[1].tolist() + pms_membs[1].tolist(), gridsize=20,
        cmap='Greens')

    plt.subplot(325)
    BV_i, V_i = np.array(field_BV_N)[msk], np.array(field_V_N)[msk]
    BV_o, V_o = np.array(field_BV_N)[~msk], np.array(field_V_N)[~msk]
    plt.scatter(BV_o, V_o, c='b', s=10, ec='k', lw=.5, alpha=.7,)
    plt.scatter(BV_i, V_i, c='orange', s=25, ec='k', lw=.5, alpha=.7,)
    plt.scatter(membs_BV, membs_V, c='r', s=30, ec='k')
    plt.gca().invert_yaxis()

    plt.subplot(326)
    plt.hexbin(
        field_BV_N + membs_BV, field_V_N + membs_V, gridsize=50,
        cmap='Greens')
    plt.gca().invert_yaxis()

    fig.tight_layout()
    plt.savefig(
        'synth_clust_out_' + str(CI) + '.png', dpi=150, bbox_inches='tight')

def CI_2(rad, cent, table_data, NN=10):

    # Select stars from cluster area (A) and field area (B)
    coord_x, coord_y = table_data['x'], table_data['y']
    coords = np.array((coord_x, coord_y)).T
    dist_cent = distance.cdist(cent, coords, 'euclidean')[0]

    # For each dimension in A, find the NN nearest neighbors in B and select the
    # one with the greatest distance
    # For each dimension in A, find the NN nearest neighbors in A and select the
    # one with the greatest distance
    data_dims = ('V', 'BV', 'pmRA', 'pmDE')
    data_arrs = []
    CI_arr = []
    for dim in data_dims:
        msk = dist_cent <= rad
        dim_a_member = table_data[dim][msk]
        dim_b_member = table_data[dim][~msk]
        dist_a_b = np.abs(dim_a_member[:, np.newaxis] - dim_b_member)
        dist_a_a = np.abs(dim_a_member[:, np.newaxis] - dim_a_member)

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
        CI_arr.append('{:.2f}'.format(CI))

        # Store for plotting
        data_arrs.append([dim_a_member, dim_b_member, CI])
    return(CI_arr, data_dims, data_arrs)

def makePlot_2(CI_coord, data_dims, data_arrs):
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
    plt.savefig("CI_analysis_" + str(CI_coord) + ".png", dpi=150, bbox_inches='tight')

if __name__ == '__main__':
    main()
