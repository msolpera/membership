
# Installation

    $ conda create -n pyupkenv numpy scikit-learn scipy astropy
    $ conda activate pyupkenv
    (pyupkenv) $ pip install kneebow





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
