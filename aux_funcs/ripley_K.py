
import numpy as np
from astropy.stats import RipleysKEstimator
from astropy.io import ascii
from astropy.table import Table
from astropy.table import Column

"""
Generate a table with the results of applying Ripley's K function over
2D uniform random fields created with different number of elements.
"""

N_max, N_runs = 2000, 1000

Kest = RipleysKEstimator(area=1)  # , x_max=1, y_max=1, x_min=0, y_min=0)
N_vals = np.arange(5, N_max)

data_dict = {}
for RK_rad in (.1, .2, .3, .4, .5, .6, .7, .8, .9):
    print(RK_rad)
    data_dict[str(RK_rad) + '_m'] = []
    data_dict[str(RK_rad) + '_s'] = []
    for N in N_vals:
        if not (N % 50):
            print(N)
        C_vals = []
        for _ in range(N_runs):
            xy_u = np.random.uniform(0., 1., (N, 2))
            C_vals.append(Kest(xy_u, (RK_rad,)))

        data_dict[str(RK_rad) + '_m'] += [np.round(np.mean(C_vals), 4)]
        data_dict[str(RK_rad) + '_s'] += [np.round(np.std(C_vals), 4)]

# Store in data file
tt = Table(data_dict)
col = Column(N_vals, name='N')
tt.add_column(col, index=0)
ascii.write(tt, "RK_data.dat", format='csv', overwrite=True)
