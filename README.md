
# Nearest neighbor clustering algorithm

- Probably related: [Nearest Neighbor Clustering: A Baseline Method for Consistent Clustering with Arbitrary Objective Functions, Bubeck & von Luxburg (2009)](http://www.jmlr.org/papers/v10/bubeck09a.html)


### kNN algorithm as a classifier
- [L8: Nearest neighbors](http://research.cs.tamu.edu/prism/lectures/pr/pr_l8.pdf)


### Voronoi diagram as a classification method
 - [Using Voronoi diagrams to improve classification performances when modeling imbalanced datasets](https://link.springer.com/article/10.1007/s00521-014-1780-0)
 - [Classification and K-nearest neighbours](https://www.inf.ed.ac.uk/teaching/courses/inf2b/learnnotes/inf2b-learn04-notes-nup.pdf)
 - [A Fuzzy Neural Network Approach to Classification Based on Proximity Characteristics of Patterns](http://www.cs.uoi.gr/~kblekas/papers/C8.pdf)



### Aux_func/generate_synth_clust
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
