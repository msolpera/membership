import mn_run
import numpy as np
from astropy.table import Table

# Read data
data = Table.read('synth_clust_out_CI_8.dat', format='ascii')

# Save data as numpy array
ID, coord_x, coord_y, pmRA, pmDE = np.array(
    [data['ID'], data['x'], data['y'], data['pmRA'], data['pmDE']])

mn_run.main(ID, coord_x, coord_y, pmRA, pmDE)
