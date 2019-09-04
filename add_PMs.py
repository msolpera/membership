
import numpy as np
from astropy.table import Table, vstack
from astropy.io import ascii
import matplotlib.pyplot as plt


def main(CI=.6):
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
    pms_membs, pms_field = generatePMs(N_membs, N_field)

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
        membs_x, membs_y, membs_V, membs_BV, x_fl, y_fl, field_BV_N, field_V_N,
        pms_membs, pms_field)

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
        cm2=.1):
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
    fm1  : float
       Minimum value for the `mean` range used for field stars.
    fm2  : float
       Maximum value for the `mean` range used for field stars.
    fstd1  : float
       Minimum value for the `stddev` range used for field stars.
    fstd2  : float
       Maximum value for the `stddev` range used for field stars.
    cm1  : float
       Minimum value for the scale used for the `stddev` of cluster stars.
    cm2  : float
       Maximum value for the scale used for the `stddev` of cluster stars.

    """

    # Generate PMs for field stars.
    mean = np.random.uniform(fm1, fm2, (2))
    # This defines the standard deviation for the field stars PMs.
    c1 = np.random.uniform(fstd1, fstd2)
    # Defined as a covariance matrix, since this is a bi-variate Gaussian.
    cov = np.eye(2) * c1
    # Generate PMs for all the field stars.
    pms_field = np.random.multivariate_normal(mean, cov, N_field)

    # Generate PMs for member stars.
    c2 = np.random.uniform(cm1, cm2)
    # Use the same covariance matrix used for the field stars, but scaled with
    # a factor 'c2'.
    cov = cov * c2

    # For the cluster members, use a `mean` taken randomly from the limits of
    # the PM values already defined for the field stars.

    # Limits in the field stars PM_RA
    range_RA = [pms_field.T[0].min(), pms_field.T[0].max()]
    # Random mean in the RA axis, within the defined range
    m1 = np.random.uniform(*range_RA)

    # Limits in the field stars PM_DEC
    range_DEC = [pms_field.T[1].min(), pms_field.T[1].max()]
    # Random mean in the DEC axis, within the defined range
    m2 = np.random.uniform(*range_DEC)

    # Mean for cluster stars
    mean_m = np.array([m1, m2])
    # Generate PMs for all the member stars.
    pms_membs = np.random.multivariate_normal(mean_m, cov, N_membs)

    return pms_membs.T, pms_field.T


def makePlot(
    membs_x, membs_y, membs_V, membs_BV, x_fl, y_fl, field_BV_N, field_V_N,
        pms_membs, pms_field):
    """
    """
    fig = plt.figure(figsize=(30, 10))

    plt.subplot(1, 3, 1)
    plt.scatter(x_fl, y_fl, c='b', s=15)
    plt.scatter(membs_x, membs_y, c='r', s=40, lw=1.5, ec='k')

    plt.subplot(1, 3, 2)
    plt.scatter(field_BV_N, field_V_N, c='b', s=15)
    plt.scatter(membs_BV, membs_V, c='r', s=40, lw=1.5, ec='k')
    plt.gca().invert_yaxis()

    plt.subplot(1, 3, 3)
    plt.scatter(pms_field[0], pms_field[1], c='b', s=15)
    plt.scatter(pms_membs[0], pms_membs[1], c='r', s=40, lw=1.5, ec='k')

    fig.tight_layout()
    plt.savefig('synth_clust_out.png', dpi=150, bbox_inches='tight')


if __name__ == '__main__':
    main()
