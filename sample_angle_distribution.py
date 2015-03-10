import numpy as np
import spherical_geo as sg
import matplotlib.pyplot as plt

num_bins = 50

pair_path = ''

pairs = np.genfromtxt(pair_path, delimiter=',', names=True)

angles = map(
             sg.calc_distance,
             zip(pairs['ra1'], pairs['dec1']),
             zip(pairs['ra2'], pairs['dec2'])
             )

counts, bin_edges = np.histogram(angles, num_bins)

cdf = np.cumsum(counts)

plt.plot(bin_edges[1:], cdf)
plt.show()