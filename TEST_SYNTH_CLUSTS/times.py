import numpy as np
import matplotlib.pyplot as plt

def main():

    UPMASK_t = 497184
    method = (
        'KNN', 'AGG', 'VOR', 'MBK', 'KMS', 'GMM', 'UPMASK\nMST', 'UPMASK')
    # Times used by the rest of the methods
    times = np.array([
        2905, 3343, 12595, 81185, 252075, 320836, 366232, UPMASK_t])
    times = times / UPMASK_t

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.set_ylabel('Time performance', fontsize=10)
    bar_plot = plt.bar(method, times, color=(
        'green', 'green', 'green', 'green', 'green', 'green', 'orange', 'red'))

    s = ('170x', '150x', '40x', '6x', '2x', '1.6x', '1.4x', '')
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)

    autolabel(bar_plot, ax, s)
    fig.tight_layout()
    plt.savefig('plots/times.png', dpi=300, bbox_inches='tight')
    # plt.figure(constrained_layout=True)
    # plt.show()


def autolabel(bar_plot, ax, s):
    for idx, rect in enumerate(bar_plot):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width() / 2., 1. * height, s[idx],
                fontsize=10, ha='center', va='bottom', rotation=0)


if __name__ == '__main__':
    main()
