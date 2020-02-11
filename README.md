
# Nearest neighbor clustering algorithm

- Probably related: [Nearest Neighbor Clustering: A Baseline Method for Consistent Clustering with Arbitrary Objective Functions, Bubeck & von Luxburg (2009)](http://www.jmlr.org/papers/v10/bubeck09a.html)


### kNN algorithm as a classifier
- [L8: Nearest neighbors](http://research.cs.tamu.edu/prism/lectures/pr/pr_l8.pdf)


### Voronoi diagram as a classification method
 - [Using Voronoi diagrams to improve classification performances when modeling imbalanced datasets](https://link.springer.com/article/10.1007/s00521-014-1780-0)
 - [Classification and K-nearest neighbours](https://www.inf.ed.ac.uk/teaching/courses/inf2b/learnnotes/inf2b-learn04-notes-nup.pdf)
 - [A Fuzzy Neural Network Approach to Classification Based on Proximity Characteristics of Patterns](http://www.cs.uoi.gr/~kblekas/papers/C8.pdf)

# memb_algor.py 
This method assumes that the stars belonging to the same stellar cluster have similar parameters, and defines the probabilities of membership of each star as a function that depends of its distance to two defined centroids.

For this, the average N-dimensional distance of each star to its n nearest neighbors is calculated. Then two percentile values for the distances, perc1 and perc2, are defined, and two groups are selected: one consisting of stars with distances less than perc1 (supposed member stars) and another with distances greater than perc2 (supposed field stars). With these stars, two N-dimensional centroids are defined as the average of their parameters,and the distance of each star from the frame to them is calculated. A function that depends on these distances is defined to give a membership probability value.

# Voronoi_v1
In this method it is assumed again that the N parameters of the member stars are similar and therefore have small Voronoi volumes associated in the N-dimensional space.

A percentile value is set, which separates the stars into two groups. Those with a Voronoi volume less than the given percentile are considered possible member stars and those with a higher volume are considered possible field stars. With each group, its Kernel Density Estimation (kd_memb, kd_field) is generated, and each star in the frame is evaluated in both. This defines two membership probability values for each star, one for the cluster (p_m) and one for the field (p_f).

The initial percentile value is modified and the corresponding KDEs are regenerated, until the maximum average value of the difference between p_m and p_f is optimized. The percentile value that optimizes this difference is chosen to generate the final probabilities.

Finally, a unique membership probability is obtained for each star by applying Bayes' theorem to the set p_m - p_f

# Voronoi_v2
The difference of this method with Voronoi_v1 is the way to select the percentile that best separates the two populations.

For each percentile value, the probability given by the Bayes Theorem is calculated, and those stars with a probability greater than 0.5 are selected. This group of stars is then compared with those of Voronoi volume less than the fixed percentile (that is, the initially assumed member stars).

Finally, that percentile value that minimizes the difference between both groups is selected and the probability found for it is considered as the final membership probability.

# Aux_func/generate_synth_clust
The purpose of this script is to generate, from a synthetic cluster, other synthetic clusters with different contamination index values.
The contamination index (CI) is a measure of the field star contamination present in the region of the stellar cluster and an CI value is defined for each stellar parameter.

CI for spatial coordinates: It is obtained as the ratio of field stars density, over the density of stars in the cluster region. 

CI for the another parameters: 
Two regions are defined. A: cluster region; B: field region.

Then, for each dimension: 
1. Select a data dimension (except coordinates)
2. For each star in A, the N nearest neighbors in B are searched in this dimension
3. Calculate the "density of neighbors" associated with that star as N / (pi * rad ** 2), where rad is the distance to the furthest neighbor
4. Obtain the "average neighbor density" of A in B in this dimension, as the average of all densities found: d_AB
5. Repeat 1, 2, 3 but now looking for the N neighbors of A in A in that same dimension.
6. Repeat 4 but with these values to obtain the "average neighbor density" of A in A in this dimension: d_AA
7. The final contamination index for this dimension is: CI = d_AB / d_AA
