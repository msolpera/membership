
# from astropy.io import ascii
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt


"""
Compare the results of a pyUPMASK run using Kmeans (25 OL runs, 25 stars per
cluster, GUMM off, KDE off) with the results from UPMASK-MST (ie: the Cantat-
Gaudin modification). This is to check wich RFR method is more adequate;
Ripley's K function or mimimum-spanning tree.'
"""

# Stars per cluster used in the Cantat-Gaudin code
stpcl = '15'  # '25'
# Folder with the results from the Cantat-Gaudin method
CG_fold = 'UPMASK_600_CG_' + stpcl
# Folder with the results from pyUPMASK using Kmeans and no GUMM or KDE
KM_fold = 'KMS_clean'

# Read metrics file for the clean Kmeans performed with pyUPMASK on 40
# synthetic clusters
pyUP_PHOT = Table.read(
    'metrics/metrics_PHOT_' + KM_fold + '_H_auto.dat', format='csv')
pyUP_PM = Table.read(
    'metrics/metrics_PM_' + KM_fold + '_H_auto.dat', format='csv')
CI_PHOT = pyUP_PHOT['CI'].data
CI_PM = pyUP_PM['CI'].data

# Read the metrics file for the 600 synthetic clusters processed with
# Cantat-Gaudin UPMASK
UPCT_PHOT = Table.read(
    'metrics/metrics_PHOT_' + CG_fold + '_H_auto.dat', format='csv')
UPCT_PM = Table.read(
    'metrics/metrics_PM_' + CG_fold + '_H_auto.dat', format='csv')

# Save as list to use the .index() method
UPCT_PHOT_names = UPCT_PHOT['Name'].tolist()
UPCT_PM_names = UPCT_PM['Name'].tolist()

metrics = ('LSR', 'BSL', 'HMS', 'MCC_5', 'TPR_5', 'PPV_5', 'MCC_9',
           'TPR_9', 'PPV_9')

m_delta_phot = []
for clst in pyUP_PHOT:
    ii = UPCT_PHOT_names.index(clst['Name'])
    # Metrics deltas.
    # Positive --> Ripley's K is better, Negative --> MST is better
    m_delta_phot.append(
        np.array([_ for _ in clst[metrics]]) -
        np.array([_ for _ in UPCT_PHOT[ii][metrics]]))
m_delta_phot = np.array(m_delta_phot)

m_delta_pm = []
for clst in pyUP_PM:
    ii = UPCT_PM_names.index(clst['Name'])
    m_delta_pm.append(
        np.array([_ for _ in clst[metrics]]) -
        np.array([_ for _ in UPCT_PM[ii][metrics]]))
m_delta_pm = np.array(m_delta_pm)

m_delta_phot_mean = np.median(m_delta_phot, 0)
m_delta_pm_mean = np.median(m_delta_pm, 0)

print("Median differences (Ripley-MST):\n")

print("PHOT")
print(("LSR={:.3f}, BSL={:.3f}, HMS={:.3f}, MCC_5={:.3f}, TPR_5={:.3f}\n" +
       "PPV_5={:.3f}, MCC_9={:.3f}, TPR_9={:.3f}, PPV_9={:.3f}").format(
           *m_delta_phot_mean))

print("\nPM")
print(("LSR={:.3f}, BSL={:.3f}, HMS={:.3f}, MCC_5={:.3f}, TPR_5={:.3f}\n" +
       "PPV_5={:.3f}, MCC_9={:.3f}, TPR_9={:.3f}, PPV_9={:.3f}").format(
           *m_delta_pm_mean))

fig, ax = plt.subplots(figsize=(15, 10))

for i in range(9):
    plt.subplot(int("33" + str(i + 1)))
    txt = r"$\Delta$PHOT={:.3f}".format(m_delta_phot_mean[i])
    txt += r", $\Delta$PM={:.3f}".format(m_delta_pm_mean[i])
    plt.title(txt, fontsize=10)
    plt.axhline(0., ls=':')
    plt.scatter(CI_PHOT, m_delta_phot[:, i], label="PHOT")
    plt.scatter(CI_PM, m_delta_pm[:, i], label="PM")
    plt.xlabel("CI")
    plt.ylabel(metrics[i])
    if i == 0:
        plt.legend()

# plt.show()
fig.tight_layout()
plt.savefig(
    'plots/RK_vs_MST_{}.png'.format(stpcl), dpi=150, bbox_inches='tight')
