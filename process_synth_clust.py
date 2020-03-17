import memb_algor
import Voronoi_v1
import Voronoi_v2
import random_memb
import P05
import UP
from aux_funcs import generate_synth_clust
import numpy as np
from astropy.table import Table
from pathlib import Path
from os import listdir
import os.path
from astropy.io import ascii

# Read files from the input folder 
def ls(ruta = Path.cwd()):
    return [arch.name for arch in Path('input').iterdir() if arch.is_file()]

def main():
    methods = ('Voronoi_v1', 'memb_algor', 'random_memb', 'P05', 'UPMASK')
    CI, CI_V, CI_BV, CI_pmRA, CI_pmDE = [], [], [], [], []
    N_m, N_f = [], []
    MI_v1, MI_ma, MI_rd, MI_p05, MI_UP = [], [], [], [], []
    p_m_v1, p_m_ma, p_m_rd, p_m_p05, p_m_UP = [], [], [], [], []
    p_f_v1, p_f_ma, p_f_rd, p_f_p05, p_f_UP = [], [], [], [], []
    arch = ls('input')
    for file_name in arch:
        CI.append(float(file_name[:3]))
        CI_V.append(float(file_name[4:8])) 
        CI_BV.append(float(file_name[9:13]))
        CI_pmRA.append(float(file_name[14:18])) 
        CI_pmDE.append(float(file_name[19:23]))
        for met in methods:
            print('Method:', met)
            print(file_name)
            if met == 'Voronoi_v1':
                ID, memb_prob = Voronoi_v1.main('input/' + file_name)
                MI, p_m, p_f, n_m, n_f = member_index(ID, memb_prob)
                MI_v1.append(MI)
                p_m_v1.append(p_m)
                p_f_v1.append(p_f)
                N_m.append(n_m)
                N_f.append(n_f)
                '''
            if met == 'Voronoi_v2':
                ID, memb_prob = Voronoi_v2.main('input/' + file_name)
                MI, p_m, p_f, n_m, n_f = member_index(ID, memb_prob)
                MI_v2.append(MI)
                p_m_v2.append(p_m)
                p_f_v2.append(p_f)
                '''
            if met == 'memb_algor':
                ID, memb_prob_ma = memb_algor.main('input/' + file_name)
                MI, p_m, p_f, n_m, n_f = member_index(ID, memb_prob_ma)
                MI_ma.append(MI)
                p_m_ma.append(p_m)
                p_f_ma.append(p_f)
            
            if met == 'random_memb':
                ID, memb_prob = random_memb.main('input/' + file_name)
                MI, p_m, p_f, n_m, n_f = member_index(ID, memb_prob)
                MI_rd.append(MI)
                p_m_rd.append(p_m)
                p_f_rd.append(p_f)
            
            if met == 'P05':
                ID, memb_prob = P05.main('input/' + file_name)
                MI, p_m, p_f, n_m, n_f = member_index(ID, memb_prob)
                MI_p05.append(MI)
                p_m_p05.append(p_m)
                p_f_p05.append(p_f)

            if met == 'UPMASK':
                ID, memb_prob = UP.main('input_up/up-' + file_name)
                MI, p_m, p_f, n_m, n_f = member_index(ID, memb_prob)
                MI_UP.append(MI)
                p_m_UP.append(p_m)
                p_f_UP.append(p_f)

            
    ascii.write([N_m, N_f, CI, CI_V, CI_BV, CI_pmRA, CI_pmDE, MI_v1, MI_ma, MI_rd, MI_p05, MI_UP], 'values.dat', names=['N_m', 'N_f', 'CI', 'CI_V', 'CI_BV', 'CI_pmRA', 'CI_pmDE', 'MI_v1', 'MI_ma', 'MI_rd', 'MI_p05', 'MI_UP'], overwrite=True)


def member_index(ID, memb_prob):
    MI = 0.
    p_m = 0.
    p_f = 0.
    n_m, n_f = 0., 0.
    for i in range(len(memb_prob)):
        if str(ID[i])[0] == '1':
            MI = MI + np.log(memb_prob[i])
            p_m = p_m + memb_prob[i]
            n_m = n_m + 1
        else:
            MI = MI + np.log(1. - memb_prob[i])
            p_f = p_f + memb_prob[i]
            n_f = n_f + 1
    p_m = p_m/n_m
    p_f = p_f/n_f
    return(MI, p_m, p_f, n_m, n_f)



'''
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
'''
if __name__ == '__main__':
    main()