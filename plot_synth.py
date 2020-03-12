import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from astropy.table import Table
import numpy as np

def main():
    data = Table.read('values.dat', format='ascii')
    data = data['N_m', 'N_f', 'CI', 'CI_V', 'CI_BV', 'CI_pmRA', 'CI_pmDE', 'MI_v1',
        'MI_v2', 'MI_ma', 'MI_rd', 'MI_p05']
    CI, MI_v1, MI_v2, MI_ma, MI_rd, MI_p05 = data['CI'], data['MI_v1'], data['MI_v2'], data['MI_ma'], data['MI_rd'], data['MI_p05']
    plt.scatter(CI-0.01, -np.log(-MI_v1), c='red', label='Voronoi_v1', s=5.)
    plt.scatter(CI-0.02, -np.log(-MI_v2), c='blue', label='Voronoi_v2', s=5.)
    plt.scatter(CI, -np.log(-MI_ma), c='green', label='Memb_algor', s=5.)
    plt.scatter(CI+0.01, -np.log(-MI_rd), c='violet', label='Random', s=5.)
    plt.scatter(CI+0.02, -np.log(-MI_p05), c='cyan', label='P05', s=5.)
    plt.xlabel('CI')
    plt.ylabel('MI')
    plt.legend()
    #plt.show()
    plt.savefig('mi_ci_out.png', dpi=150, bbox_inches='tight')

    '''
    figure(2)
    plt.scatter(CI, p_m_v1, c='red', label='Voronoi_v1', s=10.)
    plt.scatter(CI, p_m_v2, c='blue', label='Voronoi_v2', s=10.)
    plt.scatter(CI, p_m_ma, c='green', label='Memb_algor', s=20.)
    plt.scatter(CI, p_m_rd, c='violet', label='Random', s=20.)
    plt.xlabel('CI')
    plt.ylabel('Prob_memb')
    plt.legend()
    plt.savefig('prob_memb.png', dpi=150, bbox_inches='tight')
    
    figure(3)
    plt.scatter(CI, p_f_v1, c='red', label='Voronoi_v1', s=10.)
    plt.scatter(CI, p_f_v2, c='blue', label='Voronoi_v2', s=10.)
    plt.scatter(CI, p_f_ma, c='green', label='Memb_algor', s=20.)
    plt.scatter(CI, p_f_rd, c='violet', label='Random', s=20.)
    plt.xlabel('CI')
    plt.ylabel('Prob_field')
    plt.legend()
    plt.savefig('prob_field.png', dpi=150, bbox_inches='tight')
    '''
if __name__ == '__main__':
    main()