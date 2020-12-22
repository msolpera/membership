
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def main():
    n_star = ('15', '25', '50')
    lista = []
    for i, N in enumerate(n_star):
        data = ascii.read('input/metrics_UP_' + N + '.dat')
        LSR, BSL, HMS, MCC_5, TPR_5, PPV_5, MCC_9, TPR_9, PPV_9, CI = list(data['LSR']),\
            list(data['BSL']), list(data['HMS']), list(data['MCC_5']), list(data['TPR_5']), list(data['PPV_5']),\
            list(data['MCC_9']), list(data['TPR_9']), list(data['PPV_9']), np.array(data['CI'])
        list_st = (LSR + BSL + HMS + MCC_5 + TPR_5 + PPV_5 + MCC_9 + TPR_9 + PPV_9)
        lista.append(list_st)
        CI = np.repeat(CI, 9)
    lista = np.array(lista)
    # CI, st = zip(*sorted(zip(CI, list_st)))
    n15_n25 = lista[0]-lista[1]
    n15_n50 = lista[0]-lista[2]
    n25_n50 = lista[1]-lista[2]

    CI = np.array(CI)
    for l in range(0, 200, 20):
        msk = ((CI < (l+20)) & (CI >= l))
        CI[msk] = l+10
    CI = CI.astype(int)
    fig, axes = plt.subplots(1, 3, figsize=(15,5))
    aa = {'x': CI, 'y': n15_n25}
    bb = {'x': CI, 'y': n25_n50}
    cc = {'x': CI, 'y': n15_n50}
    sns.boxplot(x='x', y='y', data=aa, ax=axes[0]).set(xlabel='CI', ylabel=r'$\Delta$ metrics', title=r'$N_{15}-N_{25}$', aspect=9)
    sns.boxplot(x='x', y='y', data=bb, ax=axes[1]).set(xlabel='CI', ylabel=r'$\Delta$ metrics', title=r'$N_{25}-N_{50}$', aspect=8.7)
    sns.boxplot(x='x', y='y', data=cc, ax=axes[2]).set(xlabel='CI', ylabel=r'$\Delta$ metrics', title=r'$N_{15}-N_{50}$', aspect=7.45)
    axes[0].axhline(0, ls='--', c='red')
    axes[1].axhline(0, ls='--', c='red')
    axes[2].axhline(0, ls='--', c='red')
    fig.tight_layout()
    # plt.savefig('dmetrics_vs_ci.png', dpi=300, bbox_inches='tight')
    plt.show()


if __name__ == '__main__':
    main()
