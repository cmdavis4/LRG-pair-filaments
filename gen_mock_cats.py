import numpy as np
import math
import scipy.integrate as si

descriptor = 'centered_on_mid'
num_neighbs = 96
out_dir = '/Users/Charles/Java_Workspace/LRG_filament/'

pair_radius = 4 #Mpc
pair_z = .4
pair_id = 1
midpoint = (0,0) #RA, dec
theta = 0
delta_z_neighbs = .02
radial_bins = 4
min_rad = 3
max_rad = 12
radii = np.linspace(min_rad, max_rad, radial_bins)
neighbs_per_bin = num_neighbs / radial_bins

omega_m = .3
omega_k = 0
omega_l = .7
Dh = 3000 #h^-1 Mpc

def angDiamDistSingle(z):
    E_inv = lambda z: 1. / math.sqrt((omega_m * (1. + z)**3) + (omega_k * (1+z)**2) + omega_l)
    Dc = Dh * si.quad(E_inv, 0, z)[0]
    Dm = Dc
    Da = (Dm / (1. + z)) * (1. / .7)
    return Da

def calc_distance(ra1, dec1, ra2, dec2):
    '''Calculate the circular angular distance of two points on a sphere.'''
    lambda_diff = ra1  - ra2
    cos1 = np.cos(dec1)
    cos2 = np.cos(dec2)
    sin1 = np.sin(dec1)
    sin2 = np.sin(dec2)
    
    num = (cos2 * np.sin(lambda_diff)) ** 2.0 + (cos1 * sin2 - sin1 * cos2 * np.cos(lambda_diff)) ** 2.0
    denom = sin1 * sin2 + cos1 * cos2 * np.cos(lambda_diff)
    
    return np.arctan2(np.sqrt(num), denom)

lrg_add = angDiamDistSingle(pair_z)
lrg_ang_sep = pair_radius / lrg_add
pairs = np.zeros(12)
pairs[0] = pair_id
pairs[1] = midpoint[0]
pairs[2] = midpoint[1]
pairs[5] = lrg_ang_sep * np.cos(theta) + midpoint[0]
pairs[6] = lrg_ang_sep * np.sin(theta) + midpoint[1]
pairs[10] = lrg_ang_sep * np.cos(theta + math.pi) + midpoint[0]
pairs[11] = lrg_ang_sep * np.sin(theta + math.pi) + midpoint[1]
pairs[1:] = np.degrees(pairs[1:])
pairs[3]= pair_z
#print lrg_ang_sep
#print '%f, %f' % (pairs[5], pairs[6])

pairs = np.vstack((pairs, pairs, pairs))

neighbs = np.zeros((num_neighbs, 3))
neighbs[:,0] = pair_id
for i in range(num_neighbs):
    radius = radii[int(i / neighbs_per_bin)]
    ang_sep = radius / lrg_add
    neighbs[i, 1] = ang_sep * np.cos(2 * math.pi * (i % neighbs_per_bin) / neighbs_per_bin)
    neighbs[i, 2] = ang_sep * np.sin(2 * math.pi * (i % neighbs_per_bin) / neighbs_per_bin)
neighbs[:,1:] = np.degrees(neighbs[:,1:])

np.savetxt(out_dir + descriptor + '_pairs.cat', pairs)
np.savetxt(out_dir + descriptor + '_neighbors.cat', neighbs)
