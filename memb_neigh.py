import mn_run
import add_PMs
import numpy as np
from astropy.table import Table

for i in np.arange(0.01, 0.99, 0.01):
    CI = i
    data = add_PMs.main(CI)
    file_name = 'synth_clust_out.dat'
    mn_run.main(file_name, CI)
