import numpy as np
import spherical_geo as sg
import matplotlib.pyplot as plt
import random
import bisect as bs

num_bins = 50
num_fakes = 60000
fake_ra1_ind = 4
fake_dec1_ind = 5
fake_ra2_ind = 8
fake_dec2_ind = 9

pair_path = '/home/chadavis/catalog_creation/LRG-pair-filaments/pairs_with_pix.csv'
out_path = '/home/chadavis/catalog_creation/LRG-pair-filaments/fake_cats/sampled_cat.csv'

print 'Reading true pairs..'
pairs = np.genfromtxt(pair_path, delimiter=',', names=True)
print 'Calculating true pair CDF...'
for j in ['ra1', 'dec1', 'ra2', 'dec2']:
    pairs[j] = map(np.radians, pairs[j])

true_angles = map(
             sg.calc_distance,
             zip(pairs['ra1'], pairs['dec1']),
             [[x] for x in zip(pairs['ra2'], pairs['dec2'])]
             )

true_angles = sorted(true_angles)

p = 1. * np.arange(len(true_angles)) / (len(true_angles) - 1)

'''
plt.plot(angles, p)
plt.show()
'''

print 'Reading fake pairs...'
fake_pairs = np.loadtxt('/home/chadavis/catalog_creation/LRG-pair-filaments/fake_cats/pair_cat_LRG_Rmax24.0_rlos100.0_complete.csv', delimiter=',', skiprows=1)
#fake_pairs = fake_pairs[0:1000]
print 'Calculating fake pair CDF...'
for j in [fake_ra1_ind, fake_dec1_ind, fake_ra2_ind, fake_dec2_ind]:
    fake_pairs[:,j] = map(np.radians, fake_pairs[:,j])
fake_angles = np.array(map(
    sg.calc_distance,
    zip(fake_pairs[:,fake_ra1_ind], fake_pairs[:,fake_dec1_ind]),
    [[x] for x in zip(fake_pairs[:,fake_ra2_ind], fake_pairs[:,fake_dec2_ind])]
    ))
with_angles = np.hstack((fake_pairs, fake_angles.reshape(len(fake_angles), 1)))
with_angles = with_angles[np.argsort(with_angles[:,-1])]
sampled_pairs = []
for i in range(num_fakes):
    x = float(random.randint(0, 1000000)) / 1000000.
    if x==1.: sampled_angle = true_angles[-1]
    else: sampled_angle = true_angles[bs.bisect(p, x)][0] # + i%2?
    #print sampled_angle
    #print bs.bisect(with_angles[:,-1], sampled_angle)
    sampled_pair = with_angles[bs.bisect(with_angles[:,-1], sampled_angle)]
    sampled_pairs.append(sampled_pair)
sampled_pairs = np.vstack((sampled_pairs))
sampled_pairs = sampled_pairs[np.argsort(sampled_pairs[:,-1])]

p_fake = 1. * np.arange(len(sampled_pairs)) / (len(sampled_pairs) - 1)

print 'Saving sampled pairs...'
np.savetxt(out_path, sampled_pairs, delimiter=',')

plt.plot(true_angles, p)
plt.plot(sampled_pairs[:,-1], p_fake)
plt.show()
