
from scipy import spatial
from scipy.optimize import differential_evolution as DE
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt


def main(file_name, CI):

    # file_name = 'haf14_match.dat'  # 'synth_clust_out.dat'#'rup152_match.dat'
    print("\nProcessing: {}".format(file_name))

    # Nearest neighbors range
    n_i, n_f = 5, 25

    # Read input data
    ID, coord_x, coord_y, Vmag, BV, pmRA, pmDE = readData(file_name)

    # Normalize all input data
    stars = norm_data((coord_x, coord_y, pmRA, pmDE))

    # Process this once here, to avoid repeating multiple times inside the
    # 'for' block below.
    tree = spatial.cKDTree(stars)

    # Calculate the average distance of each star to its nearest n neighbors
    nn_avrg_dist, memb_prob, N_n = np.zeros(len(ID)), np.zeros(len(ID)), 0
    for n in range(n_i, n_f + 1, 2):

        # Average distance to the 'n' nearest neighbors
        star_n_prom = dist_neighbours(tree, n, stars)
        nn_avrg_dist += star_n_prom

        # Use the 1th percentile of the average NN distances.
        perc = np.percentile(star_n_prom, 1)

        # Select stars with average distance less than 'perc'
        select_star, select_star_dist = [], []
        for i, star in enumerate(stars):
            if star_n_prom[i] < perc:
                select_star.append(star)
                select_star_dist.append(star_n_prom[i])

        # Calculate the N-dimensional centroid of these stars as the
        # average of their positions, weighted by their associated average
        # N-N distances.
        centroid = np.average(
            select_star, axis=0, weights=1. / np.array(select_star_dist))
        print("NN: {}, centroid: ({:.3f}, {:.3f}, {:.3f}, {:.3f})".format(
            n, *centroid))

        # Calculate the probability of membership as the reciprocal of the
        # distance to the centroid
        for i, star in enumerate(stars):
            cent_dist = np.linalg.norm(star - centroid)
            memb_prob[i] += np.exp(-.5 * cent_dist)

        N_n += 1

    nn_avrg_dist, memb_prob = nn_avrg_dist / N_n, memb_prob / N_n
    nn_avrg_dist = norm_data([nn_avrg_dist], 10.)
    nn_avrg_dist = nn_avrg_dist.T[0]
    memb_prob = norm_data([memb_prob], 10.)
    memb_prob = memb_prob.T[0]

    # Using the Differential Evolution algorithm, estimate the x,y limits
    # in the upper left corner that delimitate the most probable members.
    bound_box = [[0., .5], [0.5, .999]]
    xy_lim = [.5, .5]
    crdens_frdens = boundary(
        xy_lim, nn_avrg_dist, memb_prob, coord_x, coord_y)
    for i in np.arange(0., .5, 0.005):
        for j in np.arange(0.5, .999, 0.005):
            xy_lim = [i, j]
            crdens_frdens_ij = boundary(
                xy_lim, nn_avrg_dist, memb_prob, coord_x, coord_y)
            if crdens_frdens_ij <= crdens_frdens:
                crdens_frdens = crdens_frdens_ij
                xy_lim_min = xy_lim
    xy_lim = xy_lim_min
    print("xy_lim", xy_lim)
    '''
    res = DE(
        boundary, bound_box, popsize=50, maxiter=100, disp=True,
        args=(nn_avrg_dist, memb_prob, coord_x, coord_y))
    print(res)
    xy_lim = res.x
    '''
    cent, rad, crdens, frdens = boundary(
        xy_lim, nn_avrg_dist, memb_prob, coord_x, coord_y, True)
    print(cent, rad, crdens, frdens)

    # Generate plot
    plot_memb(
        file_name, nn_avrg_dist, memb_prob, coord_x, coord_y, Vmag, BV, pmRA,
        pmDE, xy_lim, cent, rad, CI)


def readData(file_name):
    """
    Read data
    """
    data = Table.read(file_name, format='ascii')
    data = np.array([
        data['ID'], data['x'], data['y'], data['V'], data['BV'], data['pmRA'],
        data['pmDE']])
    return data


def norm_data(data_arr, nstd=3.):
    """
    Normalize input data
    """
    data_norm = []
    for arr in data_arr:
        # Mask outliers (np.nan).
        med, std = np.nanmedian(arr), np.nanstd(arr)
        dmin, dmax = med - nstd * std, med + nstd * std
        arr = np.clip(arr, dmin, dmax)

        min_array = np.nanmin(arr)
        max_array = np.nanmax(arr)
        data_norm.append((arr - min_array) / (max_array - min_array))

    return np.array(data_norm).T


def dist_neighbours(tree, n, stars):
    inx = tree.query(stars, k=n)
    NN_dist = inx[0][:, 1:]
    star_n_prom = []
    for dist in NN_dist:
        promedio = np.mean(dist)
        star_n_prom.append(promedio)
    return star_n_prom


def boundary(xyf, nn_avrg_dist, memb_prob, coord_x, coord_y, cr_flag=False):
    """
    """

    xp = np.percentile(nn_avrg_dist, 100. * xyf[0])
    yp = np.percentile(memb_prob, 100. * xyf[1])

    # Select a portion of the stars in the upper left corner
    msk = (nn_avrg_dist < xp) & (memb_prob > yp)
    xy_coord = np.array([coord_x[msk], coord_y[msk]])
    N_in = xy_coord.shape[1]

    if N_in <= 2:
        return np.inf

    # Estimate their center and radius (ie: define a "cluster region")
    cent = xy_coord.mean(1)
    rad = np.ptp(xy_coord, 1).max() * .5
    cr_area = np.pi * rad**2

    dist_2_cent = np.linalg.norm(np.array([coord_x, coord_y]).T - cent, axis=1)
    field_msk = dist_2_cent > rad
    fr_area = (np.ptp(coord_x) * np.ptp(coord_y)) - cr_area
    field_dens = field_msk.sum() / fr_area

    field_msk = dist_2_cent <= rad
    clust_reg_field_dens = max(0., field_msk.sum() - N_in) / cr_area

    if cr_flag:
        return cent, rad, clust_reg_field_dens, field_dens

    return abs(clust_reg_field_dens - field_dens)


def plot_memb(
    file_name, nn_avrg_dist, memb_prob, coord_x, coord_y, Vmag, BV, pmRA, pmDE,
        xy_lim, cent, rad, CI):
    """
    """

    memb_RA, memb_DE, memb_V, memb_BV, memb_pmRA, memb_pmDE, memb_color =\
        [[] for _ in range(7)]
    field_RA, field_DE, field_V, field_BV, field_pmRA, field_pmDE =\
        [[] for _ in range(6)]

    xp = np.percentile(nn_avrg_dist, 100 * xy_lim[0])
    yp = np.percentile(memb_prob, 100. * xy_lim[1])

    for i, x in enumerate(coord_x):
        if nn_avrg_dist[i] < xp and memb_prob[i] > yp:
            memb_RA.append(coord_x[i])
            memb_DE.append(coord_y[i])
            memb_V.append(Vmag[i])
            memb_BV.append(BV[i])
            memb_pmRA.append(pmRA[i])
            memb_pmDE.append(pmDE[i])
            memb_color.append(memb_prob[i])
        else:
            field_RA.append(coord_x[i])
            field_DE.append(coord_y[i])
            field_V.append(Vmag[i])
            field_BV.append(BV[i])
            field_pmRA.append(pmRA[i])
            field_pmDE.append(pmDE[i])

    # Define output figure size
    fig = plt.figure(figsize=(20, 20))

    plt.subplot(221)
    plt.title("N (filtered)={}".format(len(memb_RA)))
    plt.grid(ls=':', c='grey', lw=.7)
    # plt.hist(mp_norm, bins=50, density=True)
    plt.scatter(nn_avrg_dist, memb_prob, s=5)
    ymin, ymax = plt.gca().get_ylim()
    st_mean, st_std = np.median(nn_avrg_dist), np.std(nn_avrg_dist)
    xmin, xmax = max(0., min(nn_avrg_dist) - .5 * st_std),\
        st_mean + 2. * st_std
    plt.hlines(
        yp, xmin=xmin, xmax=xp, color='r', lw=2, ls='--',
        label='yp={:.2f}'.format(yp))
    plt.vlines(
        xp, ymin=yp, ymax=ymax, color='r', lw=2, ls='--',
        label='xp={:.2f}'.format(xp))
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel("Average NN distance")
    plt.ylabel('MPs')
    plt.legend()

    plt.subplot(222)
    plt.title("(x0, y0)=({:.2f}, {:.2f}) , r={:.2f}".format(*cent, rad))
    plt.grid(ls=':', c='grey', lw=.7)
    plt.scatter(memb_RA, memb_DE, s=50, c=memb_color, lw=.5,
                edgecolor='k')
    plt.colorbar(aspect=90, pad=0.01)
    plt.scatter(
        field_RA, field_DE, s=5, facecolor='none', edgecolor='k', alpha=.5)
    circle = plt.Circle(
        (cent[0], cent[1]), rad, color='r', fill=False, lw=1.5, zorder=5)
    plt.gca().add_artist(circle)
    plt.xlabel('x')
    plt.ylabel('y')

    plt.subplot(223)
    plt.grid(ls=':', c='grey', lw=.7)
    plt.scatter(memb_BV, memb_V, s=50, c=memb_color, lw=.5,
                edgecolor='k', zorder=4)
    plt.scatter(
        field_BV, field_V, marker='^', s=5, c='grey', alpha=.75, zorder=1)
    plt.xlabel('BV')
    plt.ylabel('V')
    plt.gca().invert_yaxis()

    plt.subplot(224)
    plt.grid(ls=':', c='grey', lw=.7)
    plt.scatter(memb_pmRA, memb_pmDE, s=50, c=memb_color, lw=.5,
                edgecolor='k', zorder=4)
    plt.colorbar(aspect=90, pad=0.01)
    plt.scatter(
        field_pmRA, field_pmDE, marker='^', s=5, c='grey', alpha=.75, zorder=1)
    xmean, xstd = np.mean(memb_pmRA), np.std(memb_pmRA)
    ymean, ystd = np.mean(memb_pmDE), np.std(memb_pmDE)
    plt.xlim(xmean - 3. * xstd, xmean + 3. * xstd)
    plt.ylim(ymean - 3. * ystd, ymean + 3. * ystd)
    plt.xlabel('pmRA')
    plt.ylabel('pmDE')

    fig.tight_layout()
    plt.savefig(file_name[:-4] + '_nn_' + str(CI) + '.png', dpi=150, bbox_inches='tight')


if __name__ == '__main__':
    main()
