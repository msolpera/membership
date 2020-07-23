
# Voronoi

Up until `autoperc_inner_GUMM7` all the runs use `N_membs=50`.

## autoperc_inner_GUMM
Applies the GUMM after each inner loop to reject non-members. The GUMM after the outer loop is turned off.

* UPMASK NU=25: PM 5.43 (9 wins), PHOT 0.38 (1 wins)
* UPMASK NU=50: PM 3.97 (8 wins), PHOT 1.84 (5 wins)

## autoperc_inner_GUMM2
Same as `autoperc_inner_GUMM` but adding back the GUMM after the outer loop.

* UPMASK NU=25: PM 6.04 (9 wins), PHOT 0.38 (1 wins)
* UPMASK NU=50: PM 4.51 (8 wins), PHOT 1.84 (5 wins)

## autoperc_inner_GUMM3
Same as `autoperc_inner_GUMM2` but adding a 0.05 to the GUMM probability. Marginally better than `autoperc_inner_GUMM2` in PM, much better in PHOT.

* UPMASK NU=25: PM 7.53 (9 wins), PHOT 0.43 (2 wins)
* UPMASK NU=50: PM 4.76 (8 wins), PHOT 2.27 (6 wins)

## autoperc_inner_GUMM4
Same as `autoperc_inner_GUMM2` but adding a 0.1 to the GUMM probability.

* UPMASK NU=25: PM 8.29 (9 wins), PHOT 0.49 (2 wins)
* UPMASK NU=50: PM 5.28 (8 wins), PHOT 2.11 (6 wins)

## autoperc_inner_GUMM5
Same as `autoperc_inner_GUMM2` but adding a 0.01 to the GUMM probability. Very similar to the results from `autoperc_inner_GUMM2` (expected).

* UPMASK NU=25: PM 6.20 (9 wins), PHOT 0.39 (2 wins)
* UPMASK NU=50: PM 4.52 (8 wins), PHOT 1.99 (6 wins)

## autoperc_inner_GUMM6
Same as `autoperc_inner_GUMM3` but:

1. using `rad=0.5, C_thresh=0., Kest.poisson + mode=translation`
Gives bad results

2. using `rad=0.3, C_thresh=0., Kest.poisson + mode=none`
Worse then before

3. using `rad=0.3, C_thresh=1` and values from table (not `Ktest.Poisson`)
PM results are better in summary (9 wins) but worse in combined metrics (~2.5). **Very** bad results for PHOT

4. using Dixon's test statistic Lm, with 5% critical value as C_thresh, `rad = np.array([.1, .2, .3, .4, .5])`, and `mode=none` in `Kest.Lfunction`
Similar results to `autoperc_inner_GUMM2` and `autoperc_inner_GUMM5`, slightly worse than `autoperc_inner_GUMM3`

5. using `rad = np.array([.1, .2, .3, .4, .5])` and `mode='translation` with 5% critical value
Results are far worse than before

6. using `rad = np.linspace(.01, .25, 50)` and `mode='translation` with 5% critical value
Better results than `autoperc_inner_GUMM3` in PM 5.56 (this is the only run with 9 wins in PM), slightly worse in PHOT 1.86 (6 wins)

7. using `rad = np.linspace(.01, .25, 50)` and `mode=none` with 5% critical value
Slightly worse than above (8 PM wins here)

8. using `rad = np.linspace(.01, .25, 50)` and `mode='translation` with 1% critical value (`1.68 / cl_msk.sum()`)
PM performance 6.17, PHOT performance 2.0 (similar to `autoperc_inner_GUMM5`, slightly worse than `autoperc_inner_GUMM3`) <-- **BEST**
* UPMASK NU=25: PM 7.89 (9 wins), PHOT 0.40 (2 wins)
* UPMASK NU=50: PM 6.17 (8 wins), PHOT 2.00 (6 wins)

9. using `rad = np.linspace(.001, .25, 500), mode='translation` with 1% critical value
PM performance 6.05, PHOT performance 1.98; i.e. very similar to the previous run (slightly worse)

## autoperc_inner_GUMM7
Trying to improve the Voronoi run to match UPMASK's performance when `starsPerClust_kmeans=25` is used

1. Uses `autoperc_inner_GUMM6` in mode 8., with auto `n_clusters` for Voronoi (using the `rotate()` function)
* UPMASK N=25: PM 7.23 (9 wins), PHOT 0.35 (1 win)
* UPMASK N=50: PM 6.84 (9 wins), PHOT 1.29 (4 wins)

2. Same as 1., but adding +10 to the auto `n_clusters`
* UPMASK N=25: PM 3.71 (8 wins), PHOT 0.27 (1 win)
* UPMASK N=50: PM 3.29 (7 wins), PHOT 1.25 (5 wins)

3. Uses `autoperc_inner_GUMM6` in mode 8., with `n_clusters=25`
* UPMASK N=25: PM 6.90 (9 wins), PHOT 0.37 (2 wins) <-- **BEST**
* UPMASK N=50: PM 4.45 (8 wins), PHOT 1.62 (6 wins)

4. Uses `autoperc_inner_GUMM6` in mode 8., with `n_clusters=10`
* UPMASK N=25: PM 2.82 (7 wins), PHOT 0.23 (1 win)
* UPMASK N=50: PM 2.11 (6 wins), PHOT 0.64 (3 win)

5. Uses `autoperc_inner_GUMM6` in mode 8., with `n_clusters=100`
* UPMASK N=25: PM 1.72 (6 wins), PHOT 0.33 (1 win)
* UPMASK N=50: PM 1.97 (8 wins), PHOT 1.48 (5 wins)

## voronoi_norm
Uses the configuration `autoperc_inner_GUMM6, 8.` but normalizing `clust_data` before processing with Voronoi.

1. `N_memb=50`
* UPMASK NU=25: PM 6.12 (9 wins), PHOT 0.28 (2 wins)
* UPMASK NU=50: PM 4.75 (8 wins), PHOT 1.60 (6 wins)
Doesn't give better results than the `autoperc_inner_GUMMX` runs

## agglomerative_50
Testing the configuration `autoperc_inner_GUMM6, 8.` with `AgglomerativeClustering` clustering method (`OL_runs=1, N_membs=50`).

1. linkage=ward
Old run with (L_t - rad):
* UPMASK NU=25: PM 4.30 (7 wins), PHOT 1.42 (4 wins)
* UPMASK NU=50: PM 3.11 (7 wins), PHOT 4.38 (8 wins)
New run with abs(L_t - rad):
* UPMASK NU=25: PM 4.75 (7 wins), PHOT 1.16 (4 wins)
* UPMASK NU=50: PM 3.26 (7 wins), PHOT 3.69 (7 wins)
Best for PHOT

2. linkage=single
* UPMASK NU=25: PM 1.49 (5 wins), PHOT 0.12 (0 wins)
* UPMASK NU=50: PM 1.76 (6 wins), PHOT 0.23 (1 wins)

3. linkage=average
* UPMASK NU=25: PM 5.29 (9 wins), PHOT 0.36 (1 wins)
* UPMASK NU=50: PM 4.11 (7 wins), PHOT 1.53 (5 wins)

4. linkage=complete
* UPMASK NU=25: PM 7.92 (9 wins), PHOT 0.73 (3 wins)
* UPMASK NU=50: PM 4.69 (7 wins), PHOT 2.61 (7 wins)
Best for PM

## agglomerative_25
Same as  `agglomerative_50` with `linkage=ward` and `N_membs=25`

Old run with (L_t - rad):
* UPMASK NU=25: PM 2.24 (7 wins), PHOT 1.83 (6 wins)
* UPMASK NU=50: PM 1.78 (6 wins), PHOT 3.68 (7 wins)
New run with abs(L_t - rad):
* UPMASK NU=25: PM 2.11 (7 wins), PHOT 1.49 (6 wins)
* UPMASK NU=50: PM 1.74 (6 wins), PHOT 3.36 (7 wins)
Slightly better results for PHOT under NU=25, but the `N_membs=50` run is still better.

## agglomerative_KDE
Same as  `agglomerative_50` with `clRjctMethod=kdetest, C_thresh=1`

* UPMASK NU=25: PM 20.06 (9 wins), PHOT 0.34 (1 wins)
* UPMASK NU=50: PM 10.61 (8 wins), PHOT 1.05 (5 wins)
Massively improves PM performance, at the cost of a **very poor** performance for PHOT.

## agglomerative_KDEpy
Same as  `agglomerative_KDE` with `clRjctMethod=kdetestpy`

1. `N_membs=50`
* UPMASK NU=25: PM 17.68 (9 wins), PHOT 0.29 (1 wins)
Very similar to the `agglomerative_KDE` run

2. `N_membs=25`
* UPMASK NU=25: PM 2.21 (7 wins), PHOT 0.14 (1 wins)
Worsens results considerably.

## voronoi_dist
Uses the configuration `autoperc_inner_GUMM6, 8.` but setting the `metric` in `cdist` to values other than the default `euclidean`.

1. `mahalanobis`
* UPMASK NU=25: PM 7.92 (9 wins), PHOT 0.30 (2 wins)
* UPMASK NU=50: PM 6.04 (8 wins), PHOT 1.25 (4 wins)
Bad results for PHOT.

2. `cityblock`

* UPMASK NU=25: PM 8.11 (9 wins), PHOT 0.34 (2 wins)
* UPMASK NU=50: PM 5.35 (8 wins), PHOT 1.50 (5 wins)
No improvement

3.  `chebyshev`

* UPMASK NU=25: PM 8.24 (9 wins), PHOT 0.38 (2 wins)
* UPMASK NU=50: PM 6.04 (8 wins), PHOT 2.00 (6 wins)
Almost the same as using `euclidean`

## voronoi_2idx
Uses the configuration `autoperc_inner_GUMM6, 8.` with `metric=euclidean` in `cdist`, and selecting the indexes of the first **and last** `n_clusters` elements in `idx_s`.

1. `n_clusters=auto`
* UPMASK NU=25: PM 7.23 (9 wins), PHOT 0.34 (2 wins)
* UPMASK NU=50: PM 6.42 (9 wins), PHOT 1.35 (4 wins)
Improves the PM performance for `N=50` (9 vs 8 wins), at the expense of the PHOT performance (4 vs 6 wins). Almost no improvement for `N=25`

## voronoi_kmeans
Uses the configuration `autoperc_inner_GUMM6, 8.` with `metric=euclidean` in `cdist`, `n_clusters=auto`, but using Voronoi to determine the centers and creating the clusters with `Kmeans`, I used `OL_runs=1` but I'm not sure if the results are deterministic.

* UPMASK NU=25: PM 6.17 (9 wins), PHOT 0.48 (2 wins)
* UPMASK NU=50: PM 5.53 (8 wins), PHOT 2.34 (6 wins)
No real improvement.

## kNN_10_50
Replaced the `Voronoi` algorithm for a kNN. Used 10 neighbors to estimate the density, and `N_membs=50`

* UPMASK NU=25: PM 38.89 (7 wins), PHOT -13.56 (2 wins)
Slightly better PHOT results than `autoperc_inner_GUMM6`, but slightly worst PM results.

## kNN_50_50
Used 50 neighbors and `N_membs=50`

* PM 51.11 (9 wins), PHOT -6.89 (3 wins)
Much better than `kNN_10_50` in both PHOT and PM

## kNN_25_25
Used 25 neighbors and `N_membs=25`

* PM 19.78 (7 wins), PHOT -2.22 (6 wins)
Worst results for PM but still good. Great results for PHOT: although still bellow 0, there are 6 wins.

## kNN_10_25
Used 10 neighbors and `N_membs=25`

* PM 22.00 (7 wins), PHOT -11.56 (2 wins)
Similar PM but worse PHOT than `kNN_25_25`

## kNN_50_25
Used 50 neighbors and `N_membs=25`

* PM 25.78 (7 wins), PHOT -11.56 (2 wins)
Almost equal to `kNN_10_25`

## kNN_25_50
Used 25 neighbors and `N_membs=50`

* kNN_25_50 NU=25: PM 36.22 (7 wins), PHOT -14.89 (3 wins)
Better PM, worse PHOT than `kNN_25_25` (the worst actually)

## kNN_25_auto
Used 25 neighbors and `N_membs=auto`

* PM 43.11 (7 wins), PHOT -20.22 (2 wins)
Reasonable PM, bad PHOT

## kNN_10_10
Used 10 neighbors and `N_membs=10`

* PM -8.22 (6 wins), PHOT -32.44 (2 wins)
**Horrible** performance overall

## kNN_25_mean_25
Used 25 neighbors and `N_membs=25`, but using the mean of the distances instead of the distance to the NN neighbor (as did previous runs).

* NU=25: PM 17.78 (7 wins), PHOT -2.44 (5 wins)
Very similar to `kNN_25_25`, but slightly worse PHOT

## voronoi_kde_p_25
Used Voronoi with `N_membs=25, prob_cut+=0.05` and using KDEs to assign probabilities in the [0, 1] range after the OL.

* NU=25: PM 56.00 (9 wins), PHOT -4.89 (2 wins)
Better PHOT results than the `autoperc_GUMMX` results (X=2,6) and similar PM.

## voronoi_kde_p_50
Same as `voronoi_kde_p_25` but using `N_membs=50`

* NU=25: PM 43.11 (8 wins), PHOT -8.67 (2 wins)
Worse than `voronoi_kde_p_25`

## kNN_50_50_kde
Used 50 neighbors and `N_membs=50`, with the KDE estimated probs

* NU=25: PM 49.11 (9 wins), PHOT 16.67 (7 wins)
Improves the PHOT results **a lot**

## kNN_25_25_kde
Used 25 neighbors and `N_membs=25`, with the KDE estimated probs

* NU=25: PM 40.44 (9 wins), PHOT 21.56 (8 wins)
Very good results.

## agglomerative_25_kde_p
Using `N_membs=25` with the KDE estimated probs

* NU=25: PM 54.44 (9 wins), PHOT 41.11 (9 wins)
Improves the results **a lot** compared to `agglomerative_25`

## agglomerative_50_kde_p
Using `N_membs=50` with the KDE estimated probs

* NU=25: PM 45.11 (9 wins), PHOT 20.67 (8 wins)
Improves the results for PHOT versus `agglomerative_50`. PM results are very similar, but there are more wins here.


## voronoi_newcents_50
Using Voronoi with `N_membs=50` but selecting the centers of the clusters spread across all the densities (`idxs[::step]`, where `step=N_membs`)

* NU=25: PM 50.67 (9 wins), PHOT 6.67 (5 wins)
Best results for PHOTwith Voronoi yet, and good PM results.

## voronoi_newcents_25
Same as `voronoi_newcents_50` but using Voronoi with `N_membs=25`

* NU=25: PM 44.22 (9 wins), PHOT 19.33 (8 wins)
Much better than `voronoi_newcents_50`

















# GMM

Up until `autoperc_GMM5` all the runs use the "old" RK method with `OL_runs=10, N_membs=50, RK_rad=0.5, C_thresh=2`

## autoperc_GMM
Same as `autoperc_inner_GUMM3` but using the GMM method and 10 OL runs.

* UPMASK NU=25: PM 8.09 (8 wins), PHOT 17.87 (9 wins)
* UPMASK NU=50: PM 3.77 (7 wins), PHOT 12.15 (8 wins)

## autoperc_GMM2
Same as `autoperc_GMM` but without the +0.05 added value to the GUMM prob (which means it is *similar* to `autoperc_inner_GUMM5`). <-- **BEST**

* UPMASK NU=25: PM 10.04 (9 wins), PHOT 22.67 (9 wins)
* UPMASK NU=50: PM 3.77 (7 wins), PHOT 17.61 (9 wins)

## autoperc_GMM3
Same as `autoperc_GMM` but using `covariance_type=tied` and adding +0.01 to the GUMM probability (hence equivalent to `autoperc_inner_GUMM5`)

* UPMASK NU=25: PM 16.50 (9 wins), PHOT 7.50 (8 wins)
* UPMASK NU=50: PM 5.42 (7 wins), PHOT 9.47 (7 wins)

## autoperc_GMM4
Using `covariance_type=tied` and adding +0.05 to the GUMM probability (hence equivalent to `autoperc_inner_GUMM3`)

* UPMASK NU=25: PM 17.24 (9 wins), PHOT 7.96 (8 wins)
* UPMASK NU=50: PM 5.34 (7 wins), PHOT 8.38 (7 wins)

## autoperc_GMM5
Same as `autoperc_GMM` but using 50 OL runs <-- Not sure about this (Sol)

## minibatch_50
Testing the configuration `autoperc_inner_GUMM6, 8.` with `MiniBatch` clustering method (`OL_runs=25, N_membs=50`).

* UPMASK N=25: PM 3.96 (7 wins), PHOT 2.71 (5 wins)
* UPMASK N=50: PM 2488 (7 wins), PHOT 3.17 (7 wins)
Not as good as the GMM runs, but a lot faster and still better than UPMASK: PHOT ~2 mins , PM ~45 secs (per cluster).

## minibatch_vor
`autoperc_inner_GUMM6, 8.; metric=euclidean (cdist); n_clusters=auto`, using Voronoi to determine the centers and creating the clusters with `MiniBatchKMeans`. Used `OL_runs=10`.

* UPMASK NU=25: PM 15.71 (9 wins), PHOT 1.22 (6 wins)
* UPMASK NU=50: PM 7.67 (8 wins), PHOT 4.30 (8 wins)
Much better results than `minibatch_50` for PM, and better results for PHOT.
PHOT ~6 mins , PM ~30 secs (per cluster).

























# Bad performers

## 75perc
Marks stars below the 75th percentile of the GUMM probabilities as non-members.

## autoperc
Uses the `kneebow` package to select the GUMM probability below which stars are marked as non-members.

## autoperc_5
Uses `kneebow` to select the GUMM probability, and adds 0.05 to that.

## autoperc_10
Uses `kneebow` to select the GUMM probability, and adds 0.1 to that.

## marginNmemb_autoperc
Marginalizes the N_membs parameter by running the Voronoi+GUMM method 6 times with `N_membs=(20,30,40,50,60,70)`. This requires 6 runs of the OL.

## marginC_2_autoperc
The C_thresh is varied from 0. to 3 in increments of .25.

## GUMMprobs_autoperc
Replaces all 1s with GUMM probabilities.

## norm_GUMMprobs_autoperc
Normalizes the GUMM probabilities before replacing all 1s with GUMM probabilities.

## marginC_autoperc
The C_thresh is varied from 0.5 to 2 in increments of .5.

## manualperc_1
Selects the probability as 1% (0.01)

## Nmemb20_autoperc & Nmemb70_autoperc
Same as `autoperc` but using `Nmembs=20` and `Nmembs=70`

## inner_GUMM_marginC
All these tests use the method `autoperc_inner_GUMM3`.
* `C_thresh` is varied from 0. to 3. jumping to the median C_s of the previous run.
^ Utterly fails
* `C_thresh` is varied from 1. to 10. jumping to the median C_s of the previous run.
^ Produces the **best** results for PM but the **worst** results for PHOT.
* `C_thresh` is varied from 0.5 to 2. in increases of 0.5.
^ Marginally improves results for PHOT, worsens those for PM.
* `C_thresh` is varied from 2. to 5. in increases of 0.5.
^ Same exact results as **autoperc_inner_GUMM3**
* `C_thresh` is varied from 1. to 5. in increases of 0.25.
^ The results are good for PB, really bad for PHOT
* `C_thresh` is varied from 2. to inf jumping to the median C_s of the previous run
^ Very bad results for PHOT, good results for PM

## optm_GUMM
Selects the GUMM probability cut by minimizing the difference between the (normalized) number of stars rejected and the average distance of the remaining stars to the center coordinates (obtained by the GUMM)

## voronoi_xy
Use the xy coordinates (instead of the features used to obtain the Voronoi densities) to assign the 'delta' values. Using `N_membs=50`:
* UPMASK NU=25: PM 1.53 (5 wins), PHOT 0.17 (1 wins)
* UPMASK NU=50: PM 1.75 (6 wins), PHOT 0.31 (1 wins)
Results are **terrible** because almost no "clustrer" is rejected.

## voronoi_comb
Uses the configuration `autoperc_inner_GUMM6, 8.` but combining the coordinates data with the first two dimensions of `clust_data`. This combined data is only used to obtain the Voronoi volumes, the deltas are calculated using `clust_data`.

1. `N_memb=50`
* UPMASK NU=25: PM 1.62 (7 wins), PHOT 0.20 (1 wins)
* UPMASK NU=50: PM 1.61 (7 wins), PHOT 0.80 (4 wins)
Horrible results

## optics
Testing the configuration `autoperc_inner_GUMM6, 8.` with `OPTICS` clustering method (`OL_runs=1, N_membs=50`).

* UPMASK NU=25: PM 1.26 (5 wins), PHOT 0.20 (0 wins)
* UPMASK NU=50: PM 1.34 (5 wins), PHOT 0.37 (1 wins)

## voronoi_2D
Using `autoperc_inner_GUMM6, 8., N_membs=50` but using only the first two dimensions of `clust_data` in the Voronoi algorithm. For this run I use the PM output from `autoperc_inner_GUMM6, 8.`.

* UPMASK NU=25: PM 7.89 (9 wins), PHOT 0.20 (0 wins)
* UPMASK NU=50: PM 6.17 (8 wins), PHOT 0.37 (1 wins)

## agglomerative_auto
Configuration `autoperc_inner_GUMM6, 8.` with `linkage=ward`, and Voronoi `N_membs=auto` (`OL_runs=1`).

* UPMASK NU=25: PM 6.31 (9 wins), PHOT 0.36 (1 wins)
* UPMASK NU=50: PM 5.50 (9 wins), PHOT 1.63 (5 wins)
Worsens deeply the PHOT results.

## voronoi_newcents_nauto
Auto `N_membs`

* NU=25: PM 44.22 (8 wins), PHOT -20.22 (2 wins)
Very poor results for PHOT

## voronoi_newcents_15
Same as `voronoi_newcents_50` but using Voronoi with `N_membs=15`
Not good.

## voronoi_flat
Voronoi with `N_membs=25,clRjctMethod=rkfunc` and GUMM & KDEp turned off.

* NU=25: PM -32.22 (3 wins), PHOT -52.00 (1 wins)
Horrible results