
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt


# Read synthetic cluster data
data = ascii.read('synth_clust.dat')

# Identify actual members and field stars by their IDs
membs, field = [], []
for i, st in enumerate(data['id']):
    if str(st)[0] == '1':
        membs.append(data[i])
    else:
        field.append(data[i])

aa = []
for _ in range(100):
    # Generate PMs for field stars with a large dispersion
    mean = np.random.uniform(-5., 5., (2))
    c1 = np.random.uniform(0., 5.)
    cov = np.eye(2) * c1
    pms_field = np.random.multivariate_normal(mean, cov, len(field))

    # print(pms_field.T[0].min(), pms_field.T[0].max())
    # print(pms_field.T[1].min(), pms_field.T[1].max())

    # Generate PMs for member stars with a smaller dispersion
    c2 = np.random.uniform(.1, .25)
    cov = cov * c2
    m1 = np.random.uniform(pms_field.T[0].min(), pms_field.T[0].max())
    m2 = np.random.uniform(pms_field.T[1].min(), pms_field.T[1].max())
    mean_m = np.array([m1, m2])
    mean_d = np.sqrt((
        mean[0] - mean_m[0])**2 + (mean[1] - mean_m[1])**2)
    # print(c1, c2)
    # print(mean, mean_m)
    aa.append(mean_d / c1)
    pms_membs = np.random.multivariate_normal(mean_m, cov, len(membs))

    if mean_d / c1 > 5.:
        plt.scatter(*pms_field.T, c='r')
        plt.scatter(*pms_membs.T, c='b')
        plt.show()

plt.hist(aa, 50)
plt.show()
# plt.scatter(*pms_field.T, c='r')
# plt.scatter(*pms_membs.T, c='b')
# plt.show()
