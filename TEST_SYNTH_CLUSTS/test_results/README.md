

## 75perc
Marks stars below the 75th percentile of the GUMM probabilities as non-members.

## autoperc_5
Uses `kneebow` to select the GUMM probability, and adds 0.05 to that.

## autoperc_10
Uses `kneebow` to select the GUMM probability, and adds 0.1 to that.

## autoperc
Uses the `kneebow` package to select the GUMM probability below which stars are marked as non-members.

## GUMMprobs_autoperc
Replaces all 1s with GUMM probabilities.

## marginNmemb_autoperc
Marginalizes the N_membs parameter by running the Voronoi+GUMM method 6 times with `N_membs=(20,30,40,50,60,70)`. This requires 6 runs of the OL.

## marginC_autoperc
The C_thresh is varied from 0.5 to 2 in increments of .5.

## marginC_2_autoperc
The C_thresh is varied from 0. to 3 in increments of .25.

## Nmemb20_autoperc & Nmemb70_autoperc
Same as `autoperc_auto` but using `Nmembs=20` and `Nmembs=70`

## norm_GUMMprobs_autoperc
Normalizes the GUMM probabilities before replacing all 1s with GUMM probabilities.

## manualperc_1
Selects the probability as 1% (0.01)

## autoperc_inner_GUMM
Applies the GUMM after each inner loop to reject non-members. The GUMM after the outer loop is turned off.

## autoperc_inner_GUMM2
Same as `autoperc_inner_GUMM` but adding back the GUMM after the outer loop.

## autoperc_inner_GUMM3
Same as `autoperc_inner_GUMM2` but adding a 0.05 to the GUMM probability.

## autoperc_inner_GUMM4
Same as `autoperc_inner_GUMM2` but adding a 0.1 to the GUMM probability.

## autoperc_inner_GUMM5
Same as `autoperc_inner_GUMM2` but adding a 0.01 to the GUMM probability.

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


## kneeNmemb_autoperc
Use kneebow to select the number of clusters