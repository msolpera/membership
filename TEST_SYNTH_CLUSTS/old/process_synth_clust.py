import memb_algor
from aux_funcs import generate_synth_clust
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
        data = generate_synth_clust.main(CI)
        file_name = 'synth_clust_out.dat'
        ID, coord_x, coord_y, memb_prob = memb_algor(file_name, CI)
        MI = member_index(ID, coord_x, coord_y, memb_prob)
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

def member_index(ID, coord_x, coord_y, memb_prob):
    sum_pm = 0.
    sum_memb = 0.
    sum_pf = 0.
    cent, rad = (1024, 1024), 250.
    dist_2_cent = np.linalg.norm(np.array([coord_x, coord_y]).T - cent, axis=1)
    for i in range(len(ID)):
        if (str(ID[i])[0] == '1' and dist_2_cent[i] <= rad):
            sum_pm = sum_pm + memb_prob[i]

        if (str(ID[i])[0] != '1' and dist_2_cent[i] <= rad):
            sum_pf = sum_pf + memb_prob[i]

        if str(ID[i])[0] == '1':
            sum_memb = sum_memb + 1.

    MI = (sum_pm - sum_pf)/sum_memb
    return MI