from scipy.spatial import Voronoi
from scipy.spatial import ConvexHull
from scipy import integrate
from astropy.table import Table
import numpy as np
# from sklearn.neighbors import KernelDensity
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import time
from scipy import spatial

def read_data(file_name):
    data = Table.read(file_name, format='ascii')
    data = data['ID', 'x', 'y', 'V', 'BV', 'pmRA', 'pmDE']
    # data.remove_rows(np.where([c.data for c in data.mask.itercols()])[-1])
    # msk = data['Gmag'] < 16
    # data = data[msk]
    return data['ID'], data['x'], data['y'], data['V'],\
        data['BV'], data['pmRA'], data['pmDE']

def data_norm(data_arr, nstd=3.):
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

def voronoi_volumes(points):
    v = Voronoi(points)
    vol = np.zeros(v.npoints)
    for i, reg_num in enumerate(v.point_region):
        indices = v.regions[reg_num]
        if -1 in indices: # some regions can be opened
            vol[i] = np.inf
        else:
            vol[i] = ConvexHull(v.vertices[indices]).volume
    return vol


def inx_neighbours(tree, n, stars):
    '''
    Index of the n-1 neightbours 
    '''
    dist_nn = tree.query(stars, k=n)
    inx = dist_nn[1][:, 1:]
    return inx

def main(file_name):
    # Read data
    ID, coord_x, coord_y, V, BV, pmRA, pmDE = read_data(file_name)

    # Normalize input data (without outliers)
    stars, msk_data = data_norm((coord_x, coord_y, pmRA, pmDE))
    ID = ID[msk_data]
    V = V[msk_data]
    BV = BV[msk_data]
    N_stars = stars.shape[0]
    
    # Calculate Voronoi volumes
    vol = voronoi_volumes(stars)

    # Select stars with infinite volume
    msk_v_inf = vol == np.inf
    star_vol_inf = stars[msk_v_inf]

    # Get Indexes of the 10 NN of the stars with infinite volume
    tree = spatial.cKDTree(stars[~msk_v_inf])
    inf_inx = inx_neighbours(tree, 11, star_vol_inf)

    # Assign to stars with infinite volume the average volume value of their 10 NN
    inf_vol = vol[inf_inx]
    vol[msk_v_inf] = np.median(np.ma.masked_invalid(inf_vol), axis=1)

    # Plot star coordinates with the smallest and highest Voronoi volumes (exclude stars from de frame limits)
    p_max = np.percentile(vol, 95.)
    p_min = np.percentile(vol, 2.)
    msk_field = vol >= p_max
    msk_memb = vol <= p_min
    x_field, y_field, x_memb, y_memb = stars[msk_field].T[0], stars[msk_field].T[1], stars[msk_memb].T[0], stars[msk_memb].T[1]
    print(len(x_memb))
    plt.scatter(x_memb, y_memb, s=5., label='memb')
    plt.scatter(x_field, y_field, s=5., label='field')
    plt.legend()
    plt.show()

    # Select possible member stars and field stars to generate KDE
    opt_max = np.inf
    p_dif_max = 0.
    for i in np.arange(0.2, 90, 0.1):
        p = np.percentile(vol, i) 
        msk_memb = vol <= p
        id_memb = ID[msk_memb]
        field_stars = stars[~msk_memb]
        memb_stars = stars[msk_memb]
        print(i, p, msk_memb.sum())

        # Kernel Density Estimation
        kd_field = gaussian_kde(np.array(field_stars).T)
        kd_memb =  gaussian_kde(np.array(memb_stars).T)

        # Generate probabilities
        '''
        for j in range(N_stars):
            p_memb[j] =  kd_memb.integrate_box([stars.T[0,j] - 0.05,  stars.T[1,j] - 0.05, 
                stars.T[2,j] - 0.05, stars.T[3,j] - 0.05], [stars.T[0,j] + 0.05,  stars.T[1,j] + 0.05, stars.T[2,j] + 0.05, stars.T[3,j] + 0.05])
            p_field[j] = kd_field.integrate_box([stars.T[0,j] - 0.05,  stars.T[1,j] - 0.05,
                stars.T[2,j] - 0.05, stars.T[3,j] - 0.05], [stars.T[0,j] + 0.05,  stars.T[1,j] + 0.05, stars.T[2,j] + 0.05, stars.T[3,j] + 0.05])
            end = time.time() - start
        print('int:', time.time())
        '''
        p_memb = kd_memb.evaluate(stars.T)
        p_field = kd_field.evaluate(stars.T)

        # Select the probabilities p_m, p_f that maximize the difference between them
        
        p_dif = sum(abs(p_field-p_memb))/N_stars
        if p_dif > p_dif_max:
            p_dif_max = p_dif
            index = i
            p_f = p_field
            p_m = p_memb
        print('index', index) 
        '''
        # Calculate the probability of membership by applying Bayes
        probability = 1/(1+p_field/p_memb)
        msk_prob = probability >= 0.5
        id_prob = ID[msk_prob]
        opt = len(list(set(list(id_memb)+list(id_prob))-(set(list(id_memb)) & set(list(id_prob)))))

        if opt < opt_max:
            opt_max = opt
            p_m = p_memb
            p_f = p_field
            index = i
    print(index, opt_max)
    '''
    
    probability = 1/(1+p_f/p_m)

    # Plot prob_memb vs. prob_field
    plt.scatter(p_f, p_m, s=5)
    plt.xlabel('Field probability')
    plt.ylabel('Member probability')
    plt.show()

    # Plot stars with probability greater than 0.9
    prob_memb = []
    prob_field = []
    prob_tot = 0.
    x_memb, y_memb, pmra_memb, pmde_memb, V_memb, BV_memb, prob = [], [], [], [], [], [], []
    for l, st in enumerate(stars):
        if probability[l] > 0.9:
            x_memb.append(st[0])
            y_memb.append(st[1])
            pmra_memb.append(st[2])
            pmde_memb.append(st[3])
            V_memb.append(V[l])
            BV_memb.append(BV[l])
            prob.append(probability[l])
            if str(ID[l])[0] == '1':
                prob_tot += probability[l]
            if str(ID[l])[0] != '1':
                prob_tot -= probability[l]

        if str(ID[l])[0] == '1':
            prob_memb.append(probability[l])
        if str(ID[l])[0] != '1':
            prob_field.append(probability[l])

    prob_tot = prob_tot/len(x_memb)
    print('prob_tot:', prob_tot)
    
    prob_memb_mean = np.mean(prob_memb)
    prob_field_mean = np.mean(prob_field)
    print('prob_memb_mean:', prob_memb_mean)
    print('prob_field_mean:', prob_field_mean)
    

    print(len(x_memb))
    plt.hist(probability)
    plt.show()
    plt.scatter(x_memb, y_memb, s=20., c=prob, lw=.5, edgecolor='k')
    plt.colorbar(aspect=90, pad=0.01)
    plt.savefig('coord_out.png', dpi=150, bbox_inches='tight')
    plt.show()
    plt.scatter(pmra_memb, pmde_memb, s=20., c=prob, lw=.5, edgecolor='k')
    plt.colorbar(aspect=90, pad=0.01)
    plt.savefig('pm_out.png', dpi=150, bbox_inches='tight')
    plt.show()
    plt.scatter(BV_memb, V_memb, s=20., c=prob, lw=.5, edgecolor='k')
    plt.gca().invert_yaxis()
    plt.colorbar(aspect=90, pad=0.01)
    plt.savefig('V-BV.png', dpi=150, bbox_inches='tight')
    plt.show()
       

if __name__ == '__main__':
    main('input/synth_clust_out_0.6.dat')



