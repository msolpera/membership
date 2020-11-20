import numpy as np
import matplotlib.pyplot as plt

def main():

    method = ('KNN', 'AGG', 'VOR', 'MBK', 'KMS', 'GMM', 'UPMASK')
    times = np.array([2905, 3343, 12595, 81185, 252075, 320836, 497184])
    times = times/497184
    fig, ax = plt.subplots()
    ax.set_ylabel('Time performance', fontsize=10)
    bar_plot = plt.bar(method, times, color=('green', 'green', 'green', 'green', 'green', 'green', 'red'))
    s = ('170x', '150x', '40x', '6x', '2x', '1.6x', '')
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    
    
    autolabel(bar_plot, ax, s)
    fig.tight_layout()
    plt.savefig('times.png', dpi=300, bbox_inches='tight')
    fig = plt.figure(constrained_layout=True)
    plt.show()

def autolabel(bar_plot, ax, s):
        for idx,rect in enumerate(bar_plot):
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2., 1.0*height,s[idx], fontsize=10,
                ha='center', va='bottom', rotation=0)
 

if __name__ == '__main__':
    main()