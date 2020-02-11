import memb_algor
import Voronoi_v1
import Voronoi_v2
from aux_funcs import generate_synth_clust
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from pathlib import Path
from os import listdir
import os.path
from astropy.io import ascii


# Read files from the input folder 
def ls(ruta = Path.cwd()):
    return [arch.name for arch in Path('input').iterdir() if arch.is_file()]

def main(met='Voronoi_v1'):
    arch = ls('input')
    for file_name in arch:
        print(file_name)
        if met == 'Voronoi_v1':
                ID, coord_x, coord_y, memb_prob = Voronoi_v1.main('input/' + file_name)
        if met == 'Voronoi_v2':
                ID, coord_x, coord_y, memb_prob = Voronoi_v2.main('input/' + file_name)
        if met == 'memb_algor':
                ID, coord_x, coord_y, memb_prob = memb_algor.main('input/' + file_name, CI)

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

if __name__ == '__main__':
    main()