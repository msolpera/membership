import mn_run
import add_PMs
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt


'''

'''
memb_ind = []
cont_ind = []
for i in np.arange(0.1, 0.95, 0.05):
    CI = i
    for _ in range(5):
        data = add_PMs.main(CI)
        file_name = 'synth_clust_out.dat'
        MI = mn_run.main(file_name, CI)
        memb_ind.append(MI)
        cont_ind.append(CI)

plt.figure()
plt.scatter(cont_ind, memb_ind, s=5)
plt.grid(ls=':', c='grey', lw=.7)
plt.title('Member Index vs. Contamination Index')
plt.xlabel('Cont_ind')
plt.ylabel('Memb_ind')
plt.savefig('ci_vs_mi.png', dpi=150, bbox_inches='tight')
ymin, ymax = -0.5, 1.1
xmin, xmax = 0.0, 1.0
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
