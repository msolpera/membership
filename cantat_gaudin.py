
from pathlib import Path
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import time
from pyUPMASK import readFiles
from modules.dataIO import dread, dxynorm, dwrite

# typically, runNb=100. If runNb=1 you end up with 0 or 1 for each star.
OL_runs = 25
# Threshold
T = 1
# Stars per cluster
nbStarsPerCluster = 25


def main(useCorrelations=False):
    """
    """
    out_folder = "output"
    # Create 'output' folder if it does not exist
    Path('./{}'.format(out_folder)).mkdir(parents=True, exist_ok=True)

    # Process all files inside the '/input' folder
    inputfiles = readFiles()
    for file_path in inputfiles:

        print("\n\n")
        print("Processing         : {}".format(file_path.name))
        print("Threshold          : {}".format(T))
        print("Stars per cluster  : {}".format(nbStarsPerCluster))
        start = time.time()

        # Columns names
        ID_c, x_c, y_c, data_errs = 'None', 'x', 'y', False
        if 'oc_' in file_path.name:
            # This is an UPMASK synthetic cluster
            data_cols, nbPCAcomponents =\
                ['V', 'B_V', 'U_B', 'V_I', 'J_H', 'H_K'], 4
        else:
            data_cols, nbPCAcomponents = ['pmRA', 'pmDE'], 2

        # Original data
        full_data, cl_ID, xy, data, cl_errs, data_rjct = dread(
            file_path, ID_c, x_c, y_c, data_cols, data_errs)
        # Normalize (x, y) data to [0, 1]
        xy01 = dxynorm(xy)

        probs_all = dataProcess(xy01.T, data.T, OL_runs, nbPCAcomponents)
        probs_mean = np.mean(probs_all, 0)

        # Write final data to file
        msk_data = np.array([True for _ in data])
        dwrite(
            out_folder, file_path, full_data, msk_data, probs_all, probs_mean)

        endt = time.time() - start
        print("Mins used: {:.2f}".format(endt / 60.))


def evaluate_L(fp):
    """
    Find mst of set of vertex  with Prim's algorithm.
    Total length.

    mst_edge_list = list of vertices plotable with:
                    for v in result[0]:
                            v=zip(*v)
                            plt.plot(v[0],v[1],v[2],'k-') #3D example

    m_mst_d = length of each vertex

    fp=[[0,0],[0,1],...]   list of points (list of pairs, triplets...)
    """

    def dist(somelist, i, j):
        """
        Calculates the difference between two points of a given list of points.
        """
        dim = len(somelist[0])
        s = 0
        for d in range(dim):
            s = s + (somelist[i][d] - somelist[j][d])**2
        return np.sqrt(s)

    if len(fp) == 0:
        return []
    if type(fp) is np.ndarray:
        fp = fp.tolist()

    mst_edge_list = []
    m_mst_d = []
    fp_d = np.zeros(len(fp))
    fp_pred = 1 * fp

    n = 0
    Q = fp[1:]  # Q = all points except the first one
    for p in Q:
        # ll = index of point in the original (full) list (will be 1!)
        ll = fp.index(p)
        # distance between point n (=0) and point ll of list fp
        fp_d[ll] = dist(fp, n, ll)
        fp_pred[ll] = fp[0]

    while Q != []:
        min_d = 1e20
        for p in Q:
            ll = fp.index(p)
            if fp_d[ll] < min_d:
                min_d = fp_d[ll]
                closest_fp = p
        Q.remove(closest_fp)
        m = fp.index(closest_fp)
        b0 = fp[m]
        b1 = fp_pred[m]

        m_mst_d.append(fp_d[m])
        mst_edge_list.append((b0, b1))

        for p in Q:
            ll = fp.index(p)
            d = dist(fp, m, ll)
            if d < fp_d[ll]:
                fp_d[ll] = d
                fp_pred[ll] = closest_fp
    return sum(m_mst_d)


def compute_L_uniform_squareEdge1(nbpoints):
    """
    Obtained through random realisations and fitting.
    [-0.09573605  0.81398092 -0.41282276] for std
    [-0.01479002 -0.05412924 -0.49964577] for TOTAL length
    """
    totalLengthUniform = 10**(np.polyval([
        -0.09573605, 0.81398092, -0.41282276], np.log10(nbpoints)))
    sigmaTotalLengthUniform = 10**(np.polyval(
        [-0.01479002, -0.05412924, -0.49964577], np.log10(nbpoints)))
    return totalLengthUniform, sigmaTotalLengthUniform


# def compute_L_uniform_circleRadius1(nbpoints):
#     """
#     Obtained through random realisations and fitting.
#     [-0.01735977  0.07758434 -0.03825163 -0.2310157   0.23980653
#      -0.3578292 ] for std
#     [ 0.039358   -0.32327517  1.07557236 -1.83661887  2.12329814
#      -0.51718044] for TOTAL length
#     """
#     totalLengthUniform = 10**(np.polyval([
#         0.039358, -0.32327517, 1.07557236, -1.83661887, 2.12329814,
#         -0.51718044], np.log10(nbpoints)))
#     sigmaTotalLengthUniform = 10**(np.polyval(
#         [-0.01735977, 0.07758434, -0.03825163, -0.2310157, 0.23980653,
#          -0.3578292], np.log10(nbpoints)))
#     return totalLengthUniform, sigmaTotalLengthUniform


# def readData():
    # # Read data:
    # inputFileName = 'NGC_6811_G16.xml' # sys.argv[1]
    # outputFileName = inputFileName.replace('.xml', '_withproba.csv')
    # VOTable = parse_single_table(inputFileName)

    # L = VOTable.array['l']
    # B = VOTable.array['b']
    # PAR = VOTable.array['parallax']
    # PMRA = VOTable.array['pmra']
    # PMDEC = VOTable.array['pmdec']
    # ePAR = VOTable.array['parallax_error']
    # ePMRA = VOTable.array['pmra_error']
    # ePMDEC = VOTable.array['pmdec_error']
    #
    # if useCorrelations is True:
    #     parallax_pmra_corr = VOTable.array['parallax_pmra_corr']
    #     parallax_pmdec_corr = VOTable.array['parallax_pmdec_corr']
    #     pmra_pmdec_corr = VOTable.array['pmra_pmdec_corr']

    # # organise:
    # realL = np.array(L)  # save the original l column
    # L = realL * np.cos(np.radians(np.median(B)))  # use L*cosB as coordinate
    # # NB: this only works for reasonably small fields.
    # positions = (L, B)  # used for the VETO
    # # used to make groups (can be photometric bands, or astrometric
    # # measurements, or anything really)
    # bands = (PAR, PMRA, PMDEC)

    # errors = (ePAR, ePMRA, ePMDEC)  # corresponding nominal uncertainties
    # if len(bands) != len(errors):
    #     print('You need to provide a list of associated uncertainties for' +
    #           ' each measurement.')
    #     print('Check that bands and errors have the same length.')
    #     print('EXITING.')
    #     sys.exit()

    # if len(bands) != 3 and useCorrelations is True:
    #     print('This implementation can only use correlations between' +
    #           ' pmra/pmdec/parallax if all three are used.')
    #     print('You should set the variable useCorrelations to False.')
    #     print('EXITING.')
    #     sys.exit()
    # if len(bands) < nbPCAcomponents:
    #     print('asked for %i PCA components, provided %i dimensions' % (
    #         nbPCAcomponents, len(bands)))
    #     # print('You asked for %i PCA components but only provided %i '
    #     #       'dimensions to work with.' % (nbPCAcomponents, len(bands)))
    #     print('EXITING.')
    #     sys.exit()

    # positions, bands = ()

    # return positions, bands


def dataProcess(positions, bands, runNb, nbPCAcomponents):

    # initially empty, will be filled with a new list for each run
    vetoresult = []
    for run in range(runNb):
        run = run + 1  # just for display
        print('')
        print('Run nb', run, '/', runNb)

        fuzzedBands = bands
        # # For every run we will use a slightly different dataset, picking
        # # new values for each star from its errir distribution
        # # Except for the first run where we use the nominal value.
        # if run == 1:
        #     fuzzedBands = bands
        # else:
        #     # WITHOUT CORRELATIONS
        #     #
        #     fuzzedBands = np.random.normal(loc=bands, scale=errors)
        #     #
        #     # WITH CORRELATIONS
        #     #
        #     # we need to redraw star by star...
        #     if useCorrelations is True:
        #         for i in range(len(positions[0])):
        #             pmratemp = PMRA[i]
        #             pmdectemp = PMDEC[i]
        #             partemp = PAR[i]
        #             e_pmratemp = ePMRA[i]
        #             e_pmdectemp = ePMDEC[i]
        #             e_partemp = ePAR[i]
        #             #
        #             temp_pmra_pmdec_corr = pmra_pmdec_corr[i]
        #             temp_parallax_pmra_corr = parallax_pmra_corr[i]
        #             temp_parallax_pmdec_corr = parallax_pmdec_corr[i]
        #             cov = np.zeros(shape=(3, 3))
        #             cov[0][0] = e_pmratemp**2
        #             cov[1][1] = e_pmdectemp**2
        #             cov[2][2] = e_partemp**2
        #             cov[0][1] = temp_pmra_pmdec_corr * e_pmratemp * e_pmdectemp
        #             cov[1][0] = cov[0][1]
        #             cov[0][2] = temp_parallax_pmra_corr * \
        #                 e_pmratemp * e_partemp
        #             cov[2][0] = cov[0][2]
        #             cov[1][2] = temp_parallax_pmdec_corr * \
        #                 e_pmdectemp * e_partemp
        #             cov[2][1] = cov[1][2]
        #             #
        #             fuzzedBands[0][i], fuzzedBands[1][i], fuzzedBands[2][i] =\
        #                 np.random.multivariate_normal(
        #                     [pmratemp, pmdectemp, partemp], cov)

        # Now we have "fuzzed bands" to work with.
        keepCleaning = True
        # We will keep cleaning until all stars left pass the veto, or until
        # there are no stars left.
        # =1 for all at first, progressively some 1s are replaced by 0s
        keepUsingForThisRun = np.ones(len(positions[0]))

        iteration = 0
        while keepCleaning is True:
            iteration = iteration + 1
            print('  Working with %i stars...' %
                  (len(keepUsingForThisRun[keepUsingForThisRun == 1])))

            print('  Applying PCA transformation...')
            fuzzedBandsForKeptPoints = np.array(
                [b[keepUsingForThisRun == 1] for b in fuzzedBands])
            # The above only has the original length in the first iteration,
            # then it is just a subset.
            normalisedBandsForPca = np.array(
                [np.array(b) / np.std(b) for b in fuzzedBandsForKeptPoints]).T
            # perform PCA:
            pca = PCA(n_components=nbPCAcomponents)
            pca.fit(normalisedBandsForPca)
            resOfPCA = pca.transform(normalisedBandsForPca)
            # NB: if you PCA from 3 to 3 this step performs nothing.

            print('  Applying k-means clustering...')
            # apply kmeans clustering to this transformed space
            nbOfClusters = int(
                np.floor(1. * len(resOfPCA) / nbStarsPerCluster)) + 1
            # because we want a random initialisation
            seedforkmeans = int(str(time.time())[-1])
            k_means = KMeans(n_clusters=nbOfClusters, init='random',
                             random_state=seedforkmeans, verbose=0, n_init=10)
            k_means.fit(resOfPCA)
            clusterPred = np.array(k_means.predict(resOfPCA))
            clustersYES = []
            clustersNO = []

            # # check if each one of those kmeans clusters is spatially more
            # # concentrated than a random distribution:
            # nbOfKMeansCluster = len(sorted(set(clusterPred)))

            # the following line is awful but just takes the subset of stars
            # that have not yet been rejected
            # as rows in a zipped list of the star position
            positionsOfRemainingStars = np.array(
                list(zip(*positions)))[keepUsingForThisRun == 1]
            for cnb in sorted(set(clusterPred)):
                nbInThisKMeansCluster = len(clusterPred[clusterPred == cnb])

                if nbInThisKMeansCluster > 3:

                    # observed for that group of stars:
                    totMstLength = evaluate_L(
                        positionsOfRemainingStars[clusterPred == cnb])

                    # expected for a random distribution of the same number of
                    # stars:
                    # radius = (max(positions[1]) - min(positions[1])) / 2.
                    # meanLi_1, stdLi_1 = compute_L_uniform_circleRadius1(
                    #     nbInThisKMeansCluster)  # <- for a circle of radius 1
                    # meanLi = meanLi_1 * radius
                    # stdLi = stdLi_1 * radius

                    # For a square of edge length 1
                    meanLi, stdLi = compute_L_uniform_squareEdge1(
                        nbInThisKMeansCluster)

                    # stringToPrint = (
                    #     '  Run {}.{} kmeans group {}/{} with {} stars\n'
                    #     '    length={:.3f}, Lambda= {:.3f}, random if '
                    #     '{:.3f} +- {:.3f} (T={})').format(
                    #         run, iteration, cnb + 1, nbOfKMeansCluster,
                    #         nbInThisKMeansCluster, totMstLength,
                    #         1. * (meanLi - totMstLength) / stdLi, meanLi,
                    #         T * stdLi, T)
                    # stringToPrint = "Cluster {} ({} stars)".format(
                    #     cnb + 1, nbInThisKMeansCluster)
                    if (totMstLength < meanLi - T * stdLi):
                        # print(stringToPrint + ' CLUSTERED!')
                        clustersYES.append(cnb)
                    else:
                        # print(stringToPrint + ' NOT CLUSTERED...')
                        clustersNO.append(cnb)

                else:
                    # # if too few stars, automatically consider the group not
                    # # clustered.
                    # stringToPrint = '  Run %i.%i kmeans group %2i/%i with %2i stars, not enough stars, considered NOT CLUSTERED...' % (
                    #     run, iteration, cnb + 1, nbOfKMeansCluster, nbInThisKMeansCluster)
                    # print(stringToPrint)
                    clustersNO.append(cnb)

            # The following takes advantage of the fact that even if resOfPCA
            # and clusterPred are shorter than keepUsingForThisRun
            # (because keepUsingForThisRun never changes length but the other
            # two keep getting smaller at each iteration)
            # the order is the same, keepUsingForThisRun can be directly
            # matched to clusterPred to decide which stars remain "1" and
            # which have to turn "0".
            indexAmongKept = 0
            for ii, flag in enumerate(keepUsingForThisRun):
                if flag == 1:
                    if clusterPred[indexAmongKept] in clustersNO:
                        keepUsingForThisRun[ii] = 0
                    else:
                        pass
                    indexAmongKept = indexAmongKept + 1

            # if all the stars have been discarded, or if they all passed the
            # veto, then this run is finished
            if len(clustersYES) == 0 or len(clustersYES) == nbOfClusters:
                print('  \nEND OF RUN', run)
                # store the result of which stars were kept:
                vetoresult.append(keepUsingForThisRun)
                keepCleaning = False
            else:
                print('  \nSOME GROUPS ARE CLUSTERED but not all, iterating' +
                      ' again in run', run)

    return vetoresult

    # # Now write the computed probability to the original file:
    # finalproba = sum(vetoresult) / runNb
    # finalprobastring = ['%.3f' % (foo) for foo in finalproba]
    # temptable = VOTable.to_table()  # so we can manipulate it more easily

    # # try to add a column with the requested name:
    # # try:
    # proba_col = MaskedColumn(data=finalproba, name=nameofprobacolumn)
    # # now this is a Table object, that we cannot save as XML yet.
    # temptable.add_column(proba_col)
    # VOTablethatcanbewrittenout = from_table(temptable)  # now we can!
    # # writeto(VOTablethatcanbewrittenout, outputFileName)
    # table = VOTablethatcanbewrittenout.get_first_table().to_table()
    # ascii.write(table, outputFileName, format='csv', overwrite=True)
    # # # if the name already exists then it will crash, so we replace it
    # # # instead:
    # # except:
    # #     temptable[nameofprobacolumn] = finalproba
    # #     VOTablethatcanbewrittenout = from_table(temptable)
    # #     writeto(VOTablethatcanbewrittenout, outputFileName)


if __name__ == '__main__':
    main()
