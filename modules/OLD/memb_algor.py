
import os
from scipy import spatial
# from scipy.optimize import differential_evolution as DE
from astropy.table import Table
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def main(file_name):
    """
    Calls the container functions.
    Define membership probability as the function that depends of the distances to two defined centroids.

    Parameters
    ---------------
    firts : string        
        name of the file that contains the data to read.
    second : float
        in the case of a synth clust, the value of the contamination index.

    Returns
    ---------------
    var1 : list
        star ID.
    var2 and var3 : list
        star coordinates.
    var 4 : array
        membership probability value obteined by the method for each star.  

    """

    print("\nProcessing: {}".format(file_name[6:-4]))

    # Nearest neighbors range
    n_i, n_f = 5, 25

    # Read input data
    ID, coord_x, coord_y, Vmag, BV, pmRA, pmDE = readData(file_name)

    # Normalize all input data
    stars, msk_data = norm_data((pmRA, pmDE))
    N_stars = stars.shape[0]

    # Process this once here, to avoid repeating multiple times inside the
    # 'for' block below.
    tree = spatial.cKDTree(stars)

    # Calculate the average distance of each star to its nearest n neighbors
    nn_avrg_dist, memb_prob, N_n, d1, d2 = np.zeros(N_stars), np.zeros(N_stars), 0, np.zeros(N_stars), np.zeros(N_stars)
    cent_dist_1, cent_dist_2 = np.zeros(N_stars), np.zeros(N_stars)
    for n in range(n_i, n_f + 1, 2):

        # Average distance to the 'n' nearest neighbors
        star_n_prom = dist_neighbours(tree, n, stars)
        nn_avrg_dist += star_n_prom

        # Use the Xth percentile of the average NN distances.
        perc1, p = 0., 1.
        p2 = 99.
        while perc1 <= 0.:
            perc1 = np.percentile(star_n_prom, p)
            p += 1.
        perc2 = np.percentile(star_n_prom, p2)

        # Select stars with average distance less than 'perc1'
        select_star_1, select_star_dist_1 = [], []
        for i, star in enumerate(stars):
            if star_n_prom[i] <= perc1:
                select_star_1.append(star)
                select_star_dist_1.append(star_n_prom[i])
        
        # Select stars with average distance greather than 'perc2'
        select_star_2, select_star_dist_2 = [], []
        for i, star in enumerate(stars):
            if star_n_prom[i] >= perc2:
                select_star_2.append(star)
                select_star_dist_2.append(star_n_prom[i])

        # Calculate the two N-dimensional centroid of these stars as the
        # average of their positions, weighted by their associated average
        # N-N distances.

        w1 = np.exp(-.5 * np.array(select_star_dist_1))
        centroid_1 = np.average(select_star_1, axis=0, weights=w1)
        print("NN: {}, centroid: {}".format(n, centroid_1))

        w2 = np.exp(-.5 * np.array(select_star_dist_2))
        centroid_2 = np.average(select_star_2, axis=0, weights=w2)
        print("NN: {}, centroid: {}".format(n, centroid_2))

        # Calculate the probability of membership as the function below
        
        for i in range(len(stars)):
            cent_dist_1[i] += np.linalg.norm(stars[i] - centroid_1)
            cent_dist_2[i] += np.linalg.norm(stars[i] - centroid_2)
            memb_prob[i] += np.exp(-.5 * cent_dist_1[i]/cent_dist_2[i])
        N_n += 1

    d1, d2 = cent_dist_1.T, cent_dist_2.T
    data_n, _ = norm_data([memb_prob, nn_avrg_dist], 100.)
    memb_prob, nn_avrg_dist = data_n.T
    for i in range(len(stars)):
        if memb_prob[i] <= 0.01:
            memb_prob[i] = 0.01
        if memb_prob[i] == 1:
            memb_prob[i] = 0.99
    

    # Estimate the x,y limits in the upper left corner that 
    # delimitate the most probable members.
    xy_lim = [.5, .5]
    crdens_frdens = boundary(
        xy_lim, nn_avrg_dist, memb_prob, coord_x[msk_data], coord_y[msk_data])
    for i in np.arange(0., .5, 0.005):
        for j in np.arange(0.5, .999, 0.005):
            xy_lim = [i, j]
            crdens_frdens_ij = boundary(
                xy_lim, nn_avrg_dist, memb_prob, coord_x[msk_data],
                coord_y[msk_data])
            if crdens_frdens_ij <= crdens_frdens:
                crdens_frdens = crdens_frdens_ij
                xy_lim_min = xy_lim
    xy_lim = xy_lim_min
    print("xy_lim", xy_lim)

    # Split cluster members and field stars populations.
    memb_ID, memb_x, memb_y, memb_V, memb_BV, memb_pmRA, memb_pmDE,\
        memb_color, field_x, field_y, field_V, field_BV, field_pmRA,\
        field_pmDE = membSplit(
            ID, coord_x, coord_y, Vmag, BV, pmRA, pmDE, msk_data,
            nn_avrg_dist, memb_prob, xy_lim)

    cent, rad, crdens, frdens = boundary(
        xy_lim, nn_avrg_dist, memb_prob, coord_x[msk_data], coord_y[msk_data],
        True)
    print(cent, rad, crdens, frdens)

    # Store the selected member stars in an output file.
    # storeMembs(file_name, memb_ID, CI)

    # Generate plot
    #plot_memb(
    #    file_name, CI, memb_prob, nn_avrg_dist, memb_x, memb_y, memb_V,
    #    memb_BV, memb_pmRA, memb_pmDE, memb_color, field_x, field_y, field_V,
    #    field_BV, field_pmRA, field_pmDE, xy_lim, cent, rad, d1, d2)
    
    return ID, memb_prob

'''
def readData(file_name):
    """
    Read data
    """
    data = Table.read(file_name, format='ascii')
    data = np.array([
        data['DR2Name'], data['_x'], data['_y'], data['Gmag'], data['BP-RP'], data['pmRA'],
        data['pmDE']])
    return data
'''
def readData(file_name):
    """
    Read data

    Parameter
    ---------------
    String : 
        name of the file that contains the data to read.
    
    Returns
    --------------
    var1 : list
        star id.
    var2 and var3 : list
        star coordinates 
    var4 and var5 : list
        magnitudes and colour indices
    var6 and var7 : list
        proper motions
    """
def readData(file_name):
    data = Table.read(file_name, format='ascii')
    data = data['ID', 'x', 'y', 'V', 'BV', 'pmRA', 'pmDE']
    # data.remove_rows(np.where([c.data for c in data.mask.itercols()])[-1])
    # msk = data['Gmag'] < 16
    # data = data[msk]
    return data['ID'], data['x'], data['y'], data['V'], data['BV'], data['pmRA'], data['pmDE']


def norm_data(data_arr, nstd=3.):
    """
    Normalize input data

    Parameters
    --------------
    first : list
        data to normalize.
    second : float
        maximum standard deviation considered

    Returns
    --------------
    var1 : array
        normalized data
    var2 :
        masked data

    """
    msk_all = []
    for arr in data_arr:
        # Mask outliers (np.nan).
        med, std = np.nanmedian(arr), np.nanstd(arr)
        dmin, dmax = med - nstd * std, med + nstd * std
        msk = (arr > dmin) & (arr < dmax)
        msk_all.append(msk.data)

    msk_data = np.logical_and.reduce(msk_all) 

    data_norm = []
    for arr in data_arr:
        min_array = np.nanmin(arr[msk_data])
        max_array = np.nanmax(arr[msk_data])
        data_norm.append((arr[msk_data] - min_array) / (max_array - min_array))

    return np.array(data_norm).T, msk_data


def dist_neighbours(tree, n, stars):
    """
    Calculate the distance from each star to its nearest 'n' neighbors.

    Parameters:
    ---------------
    first :

    second : int
        number of nearest neighbours to take
    third : 
        normalized data
    """
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

    # xp = np.percentile(nn_avrg_dist, 100. * xyf[0])
    # yp = np.percentile(memb_prob, 100. * xyf[1])
    xp, yp = xyf

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


def membSplit(
    ID, coord_x, coord_y, Vmag, BV, pmRA, pmDE, msk_data, nn_avrg_dist,
        memb_prob, xy_lim):
    """
    """
    msk_ID, msk_RA, msk_DE, msk_V, msk_BV, msk_pmRA, msk_pmDE,\
        msk_Plx = [[] for _ in range(8)]
    field_RA, field_DE, field_V, field_BV, field_pmRA, field_pmDE, field_Plx =\
        [[] for _ in range(7)]

    for i, flag in enumerate(msk_data):
        if bool(flag) is True:
            msk_ID.append(ID[i])
            msk_RA.append(coord_x[i])
            msk_DE.append(coord_y[i])
            msk_V.append(Vmag[i])
            msk_BV.append(BV[i])
            msk_pmRA.append(pmRA[i])
            msk_pmDE.append(pmDE[i])
            # msk_Plx.append(Plx[i])
        else:
            field_RA.append(coord_x[i])
            field_DE.append(coord_y[i])
            field_V.append(Vmag[i])
            field_BV.append(BV[i])
            field_pmRA.append(pmRA[i])
            field_pmDE.append(pmDE[i])
            # field_Plx.append(Plx[i])

    memb_ID, memb_RA, memb_DE, memb_V, memb_BV, memb_pmRA, memb_pmDE,\
        memb_Plx, memb_color = [[] for _ in range(9)]

    for i, _id in enumerate(msk_ID):
        if nn_avrg_dist[i] < xy_lim[0] and memb_prob[i] > xy_lim[1]:
            memb_ID.append(_id)
            memb_RA.append(msk_RA[i])
            memb_DE.append(msk_DE[i])
            memb_V.append(msk_V[i])
            memb_BV.append(msk_BV[i])
            memb_pmRA.append(msk_pmRA[i])
            memb_pmDE.append(msk_pmDE[i])
            # memb_Plx.append(msk_Plx[i])
            memb_color.append(memb_prob[i])
        else:
            field_RA.append(msk_RA[i])
            field_DE.append(msk_DE[i])
            field_V.append(msk_V[i])
            field_BV.append(msk_BV[i])
            field_pmRA.append(msk_pmRA[i])
            field_pmDE.append(msk_pmDE[i])
            # field_Plx.append(msk_Plx[i])

    return memb_ID, memb_RA, memb_DE, memb_V, memb_BV, memb_pmRA, memb_pmDE,\
        memb_color, field_RA, field_DE, field_V, field_BV, field_pmRA,\
        field_pmDE


def storeMembs(file_name, memb_ID, CI):
    """
    """
    # If the 'output/' folder does not exist, create one.
    if not os.path.exists('output'):
        os.makedirs('output')

    out_name = 'output/' + file_name[6:-4] + '_nn_' + str(CI) + '.dat'
    ascii.write(
        [memb_ID], out_name, names=['ID'], overwrite=True)


def plot_memb(
    file_name, CI, memb_prob, nn_avrg_dist, memb_RA, memb_DE, memb_V, memb_BV,
    memb_pmRA, memb_pmDE, memb_color, field_RA, field_DE, field_V,
        field_BV, field_pmRA, field_pmDE, xy_lim, cent, rad, d1, d2):
    """
    """
    # Define output figure size
    fig = plt.figure(figsize=(20, 50))
    gs = gridspec.GridSpec(10, 4)

    ax = plt.subplot(gs[0:2, 0:2])
    ax.minorticks_on()
    plt.title("N (filtered)={}".format(len(memb_RA)))
    plt.grid(ls=':', c='grey', lw=.7)
    # plt.hist(mp_norm, bins=50, density=True)
    plt.scatter(nn_avrg_dist, memb_prob, s=5)
    ymin, ymax = plt.gca().get_ylim()
    st_mean, st_std = np.median(nn_avrg_dist), np.std(nn_avrg_dist)
    xmin, xmax = max(0., min(nn_avrg_dist) - .5 * st_std),\
        st_mean + 2. * st_std
    xp = np.percentile(nn_avrg_dist, 100 * xy_lim[0])
    yp1 = np.percentile(memb_prob, 100. * xy_lim[1])
    # yp2 = np.percentile(memb_prob, 100. * xy_lim[2])
    plt.hlines(
        xy_lim[1], xmin=xmin, xmax=xy_lim[0], color='r', lw=2, ls='--',
        label='MP={:.2f} (yp1={:.2f})'.format(xy_lim[1], yp1))
    # plt.hlines(
    #     xy_lim[2], xmin=xmin, xmax=xy_lim[0], color='r', lw=2, ls='--',
    #     label='MP={:.2f} (yp2={:.2f})'.format(xy_lim[2], yp2))
    plt.vlines(
        xy_lim[0], ymin=xy_lim[1], ymax=1., color='r', lw=2, ls='--',
        label='NN dist={:.2f} (xp={:.2f})'.format(xy_lim[0], xp))
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel("Average NN distance")
    plt.ylabel('MPs')
    plt.legend()

    plt.subplot(gs[0:2, 2:4])
    plt.title("distance to centroid 2 vs centroid 1")
    plt.grid(ls=':', c='grey', lw=.7)
    plt.scatter(d2, d1, s=5)
    plt.xlabel('Distance to centroid of great distances')
    plt.ylabel('Distance to centroid of small distances')

    plt.subplot(gs[2:4, 0:2])
    plt.title("(x0, y0)=({:.2f}, {:.2f}) , r={:.2f}".format(*cent, rad))
    plt.grid(ls=':', c='grey', lw=.7)
    plt.scatter(memb_RA, memb_DE, s=50, c=memb_color, lw=.5, edgecolor='k')
    plt.colorbar(aspect=90, pad=0.01)
    plt.scatter(
        field_RA, field_DE, s=5, facecolor='none', edgecolor='k', alpha=.5)
    circle = plt.Circle(
        (cent[0], cent[1]), rad, color='r', fill=False, lw=1.5, zorder=5)
    plt.gca().add_artist(circle)
    plt.xlabel('x')
    plt.ylabel('y')

    plt.subplot(gs[2:4, 2:4])
    # plt.title("(x0, y0)=({:.2f}, {:.2f}) , r={:.2f}".format(*cent, rad))
    plt.grid(ls=':', c='grey', lw=.7)
    plt.scatter(memb_RA, memb_DE, s=50, c=memb_color, lw=.5, edgecolor='k')
    plt.colorbar(aspect=90, pad=0.01)
    plt.xlabel('x')
    plt.ylabel('y')

    # IN PLACE FOR THE PARALAX
    # plt.subplot(gs[4:6, 0:2])
    # plt.grid(ls=':', c='grey', lw=.7)
    # plt.scatter(memb_Plx, memb_V, s=50, c=memb_color, lw=.5,
    #             edgecolor='k', zorder=4)
    # plt.colorbar(aspect=90, pad=0.01)
    # plt.scatter(
    #     field_Plx, field_V, marker='^', s=5, c='grey', alpha=.75, zorder=1)
    # xmean, xstd = np.mean(field_Plx), np.std(field_Plx)
    # plt.xlim(max(-1., xmean - 3. * xstd), xmean + 3. * xstd)
    # plt.gca().invert_yaxis()
    # plt.xlabel('Plx')
    # plt.ylabel('Mag')

    # plt.subplot(gs[4:6, 2:4])
    # plt.grid(ls=':', c='grey', lw=.7)
    # plt.scatter(
    #     field_Plx, field_V, marker='^', s=5, c='grey', zorder=1)
    # plt.scatter(memb_Plx, memb_V, s=50, c=memb_color, lw=.5, alpha=.75,
    #             edgecolor='k', zorder=4)
    # plt.colorbar(aspect=90, pad=0.01)
    # xmean, xstd = np.mean(memb_Plx), np.std(memb_Plx)
    # plt.xlim(xmean - 3. * xstd, xmean + 3. * xstd)
    # plt.gca().invert_yaxis()
    # plt.xlabel('Plx')
    # plt.ylabel('Mag')

    plt.subplot(gs[6:8, 0:2])
    plt.grid(ls=':', c='grey', lw=.7)
    plt.scatter(memb_pmRA, memb_pmDE, s=50, c=memb_color, lw=.5,
                edgecolor='k', zorder=4)
    plt.colorbar(aspect=90, pad=0.01)
    plt.scatter(
        field_pmRA, field_pmDE, marker='^', s=5, c='grey', alpha=.75, zorder=1)
    xmean, xstd = np.mean(field_pmRA), np.std(field_pmRA)
    ymean, ystd = np.mean(field_pmDE), np.std(field_pmDE)
    plt.xlim(xmean - 3. * xstd, xmean + 3. * xstd)
    plt.ylim(ymean - 3. * ystd, ymean + 3. * ystd)
    plt.xlabel('pmRA')
    plt.ylabel('pmDE')

    plt.subplot(gs[6:8, 2:4])
    plt.grid(ls=':', c='grey', lw=.7)
    plt.scatter(
        field_pmRA, field_pmDE, marker='^', s=5, c='grey', zorder=1)
    plt.scatter(memb_pmRA, memb_pmDE, s=50, c=memb_color, lw=.5, alpha=.75,
                edgecolor='k', zorder=4)
    plt.colorbar(aspect=90, pad=0.01)
    xmean, xstd = np.mean(memb_pmRA), np.std(memb_pmRA)
    ymean, ystd = np.mean(memb_pmDE), np.std(memb_pmDE)
    plt.xlim(xmean - 3. * xstd, xmean + 3. * xstd)
    plt.ylim(ymean - 3. * ystd, ymean + 3. * ystd)
    plt.xlabel('pmRA')
    plt.ylabel('pmDE')

    plt.subplot(gs[8:10, 0:2])
    plt.grid(ls=':', c='grey', lw=.7)
    plt.scatter(memb_BV, memb_V, s=50, c=memb_color, lw=.5,
                 edgecolor='k', zorder=4)
    plt.scatter(
        field_BV, field_V, marker='^', s=5, c='grey', alpha=.75, zorder=1)
    plt.xlabel('BV')
    plt.ylabel('V')
    plt.gca().invert_yaxis()

    plt.subplot(gs[8:10, 2:4])
    plt.grid(ls=':', c='grey', lw=.7)
    plt.scatter(memb_BV, memb_V, s=50, c=memb_color, lw=.5,
                edgecolor='k', zorder=4)
    # plt.scatter(
    #     field_BV, field_V, marker='^', s=5, c='grey', alpha=.75, zorder=1)
    plt.xlabel('BV')
    plt.ylabel('V')
    plt.gca().invert_yaxis()

    fig.tight_layout()

    # Output image name
    out_name = 'output/' + file_name[6:-4] + '_nn_' + str(CI) + '.png'
    plt.savefig(out_name, dpi=150, bbox_inches='tight')


if __name__ == '__main__':
    main()
