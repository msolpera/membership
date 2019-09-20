import mn_run
import numpy as np
from astropy.table import Table

# Read data
data = Table.read('haf14_match.dat', format='ascii')

# Save data as numpy array
ID, coord_x, coord_y, RA, DE, pmRA, pmDE = np.array(
    [data['id'], data['x'], data['y'], data['RA_ICRS'], data['DE_ICRS'],
        data['pmRA'], data['pmDE']])

mn_run.main(ID, RA, DE, pmRA, pmDE)