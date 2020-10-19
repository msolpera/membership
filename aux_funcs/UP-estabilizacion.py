from pathlib import Path
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt


def main():
    inputfiles = readFiles('input')

    stats = {'LSR': [], 'BSL': [], 'HMS': [], 'MCC_5': [], 'TPR_5': [],
             'PPV_5': [], 'MCC_9': [], 'TPR_9': [], 'PPV_9': []}
    rng = 0.025

    # For each synthetic cluster
    for arch in inputfiles:
        print(arch)
        # Read data from synthetic cluster's file
        data = ascii.read('input/' + arch)

        for sta in stats.keys():
            # Data for this statistic
            s = data[sta]
            break_flag = False
            for i, sn in enumerate(s):
                # Break before reaching the end of the column
                if i == len(s) - 5:
                    break
                # If the values are within range for 5 continuous steps,
                # store how many runs it took this statistic to achieve
                # convergence
                _min, _max = sn - rng, sn + rng
                if (s[i + 1] <= _max) & (s[i + 1] >= _min):
                    if (s[i + 2] <= _max) & (s[i + 2] >= _min):
                        if (s[i + 3] <= _max) & (s[i + 3] >= _min):
                            if (s[i + 4] <= _max) & (s[i + 4] >= _min):
                                if (s[i + 5] <= _max) & (s[i + 5] >= _min):
                                    stats[sta].append(data['n_run'][i + 5])
                                    break_flag = True
                                    break

            # If no convergence was found for this file, assign a large value
            if break_flag is False:
                stats[sta].append(40)

    fig, ax = plt.subplots(figsize=(6, 6))
    # plot the cumulative histogram
    for sta, lst in stats.items():
        nbins = np.arange(6, 50, 2)
        n, bins, patches = ax.hist(
            lst, bins=nbins, density=True, histtype='step', cumulative=True,
            label=sta)
    #     bins = bins[:-1]
    #     msk = n >= 0.95
    #     bins_msk = bins[msk]
    #     n_break[i] = bins_msk[0]
    #     ax.plot(n_break[i], 0.95, marker="x", color="red")
    # n_mean = np.mean(n_break)
    # n_max = np.max(n_break)

    ax.axhline(y=0.85, color='r', linestyle='--')
    # ax.grid(True)
    ax.legend(loc='lower right')
    # ax.set_title('n_mean = {:.1f}, n_max = {}'.format(n_mean, n_max), fontsize=10)
    ax.set_xlabel('Outer loop runs', fontsize=10)
    ax.set_ylabel('Convergence percentil', fontsize=10)
    x_int = [5, 10, 15, 20, 25]
    plt.xticks(x_int, fontsize=10)
    plt.yticks(fontsize=10)
    plt.xlim(5, 25)

    fig.tight_layout()
    plt.savefig('pyUP_MiniB010_convergence.png', dpi=300, bbox_inches='tight')
    # fig = plt.figure(constrained_layout=True)
    # plt.show()


def readFiles(name):
    """
    Read files from the input folder
    """
    return [arch.name for arch in Path(name).iterdir() if arch.is_file()]

if __name__ == '__main__':
    main()