
import numpy as np
from astropy.table import Table, vstack
from astropy.io import ascii
import matplotlib.pyplot as plt


def main(CI):
    """

    xy_range : float
      The x,y range of the total frame in pixels.
    cl_area: float
      The area of the cluster region. Set to be a circle of r=250 px.
    CI : float < 1.
      The contamination index. Set to be 0.6 to match Haffner 14.
    """

    # Hardcoded parameters: cluster radius and frame lengths in x & y.
    cl_rad, xy_range = 250., 2048.

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
    ascii.write(
        vstack([tabl_cl, tabl_fl]), 'synth_clust_out.dat', overwrite=True)

    # Generate plot
    makePlot(
        CI, membs_x, membs_y, membs_V, membs_BV, x_fl, y_fl, field_BV_N,
        field_V_N, pms_membs, pms_field, dmean)

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
    N_field_in_clreg = N_membs / (1. / CI - 1.)

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
        cm2=.1, nstd=2.):
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
    CI, membs_x, membs_y, membs_V, membs_BV, x_fl, y_fl, field_BV_N, field_V_N,
        pms_membs, pms_field, dmean):
    """
    """
    fig = plt.figure(figsize=(10, 15))
    plt.suptitle("CI={:.2f}, d_mean={:.2f}".format(CI, dmean), y=1.02)

    plt.subplot(321)
    plt.scatter(x_fl, y_fl, c='b', ec='k', lw=.5, alpha=.7, s=10)
    plt.scatter(membs_x, membs_y, c='r', s=40, ec='k')

    plt.subplot(322)
    plt.hexbin(
        x_fl.tolist() + membs_x, y_fl.tolist() + membs_y, gridsize=20,
        cmap='Greens')

    plt.subplot(323)
    plt.scatter(
        pms_field[0], pms_field[1], c='b', s=10, ec='k', lw=.5,
        alpha=.7,)
    plt.scatter(pms_membs[0], pms_membs[1], c='r', s=30, ec='k')

    plt.subplot(324)
    plt.hexbin(
        pms_field[0].tolist() + pms_membs[0].tolist(),
        pms_field[1].tolist() + pms_membs[1].tolist(), gridsize=20,
        cmap='Greens')

    plt.subplot(325)
    plt.scatter(
        field_BV_N, field_V_N, c='b', ec='k', lw=.5, alpha=.7, s=10)
    plt.scatter(membs_BV, membs_V, c='r', s=30, ec='k')
    plt.gca().invert_yaxis()

    plt.subplot(326)
    plt.hexbin(
        field_BV_N + membs_BV, field_V_N + membs_V, gridsize=50,
        cmap='Greens')
    plt.gca().invert_yaxis()

    fig.tight_layout()
    plt.savefig('synth_clust_out_' + str(CI) + '.png', dpi=150, bbox_inches='tight')


if __name__ == '__main__':
    main()
