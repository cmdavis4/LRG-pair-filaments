import numpy as np
import matplotlib.pyplot as plt
import bisect as bs
import math
import scipy.integrate as si
import numpy.linalg as LA

# local imports
from cosmology_settings import cosmology_params as cosmo_params
from cosmology_functions import angDiamDistSingle, angDiamDist
from spherical_geo import to_sphere, to_cartesian, rotate, calc_distance
from user_settings import project_path

pair_range = range(0, 9)
xbins = 20
ybins = 20
xinc = 5. / xbins
yinc = 5. / ybins
xbin_ledge = np.linspace(-5, 5, xbins + 1)
ybin_ledge = np.linspace(-5, 5, ybins + 1)
xbin_ledge = xbin_ledge[1:]
ybin_ledge = ybin_ledge[1:]
#print xbin_ledge

neighbor_id_ind = -1
neighbor_ra_ind = 11
neighbor_dec_ind = 12
neighbor_z_ind = 13
lrg_id_ind = 0
mid_ra_ind = 1
mid_dec_ind = 2
lrg_z_ind = 3
lrg1_ra_ind = 5
lrg1_dec_ind = 6
lrg2_ra_ind = 10
lrg2_dec_ind = 11

out_dir = ''
neighbor_path = '/home/chadavis/catalog_creation/LRG-pair-filaments/matches/matches_3.8.15_complete_nside16.csv'
pair_path = '/home/chadavis/catalog_creation/LRG-pair-filaments/pairs_with_pix.csv'
#neighbor_path = '/home/chadavis/catalog_creation/astro_image_processing/LRG/data/ra_dec_z.cat'
#pair_path = '/data2/scratch/pairs_6_10.txt'

print('Reading neighbors...')
neighbors = np.loadtxt(neighbor_path, delimiter=',')
ids = np.array(neighbors[:,neighbor_id_ind], dtype = np.int64)        

print('Reading pairs...')
pairs = np.loadtxt(pair_path)


objects = np.zeros(shape=(len(neighbors), 4))
objects[:,0] = neighbors[:,neighbor_id_ind] #id
objects[:,1] = neighbors[:,neighbor_ra_ind] #ra
objects[:,2] = neighbors[:,neighbor_dec_ind] #dec
objects[:,3] = neighbors[:,neighbor_z_ind] #z

print('Sorting neighbors...')
objects = np.array(sorted(objects, key = lambda x: x[0]), dtype=np.float64)
ids = sorted(ids)

xs = np.array([])
ys = np.array([])
offset = 0
grid = np.zeros(shape=(xbins, ybins))
outOfGrid = 0
total_neighbs = 0
oob_lrgs = 0

print('Creating match dictionary...')
idDict = {}
for j in range(len(objects)):
    idDict.setdefault(ids[j], []).append(j)

for i in pair_range:
    curr = pairs[i]
    lrg_id = curr[lrg_id_ind]
    ra_mid = np.radians(curr[mid_ra_ind])
    dec_mid = np.radians(curr[mid_dec_ind])
    ra_1 = np.radians(curr[lrg1_ra_ind])
    dec_1 = np.radians(curr[lrg1_dec_ind])
    ra_2 = np.radians(curr[lrg2_ra_ind])
    dec_2 = np.radians(curr[lrg2_dec_ind])
    lrg_add = angDiamDistSingle(curr[lrg_z_ind])
    
    left_ind = bs.bisect_left(ids, np.int64(lrg_id))
    right_ind = bs.bisect(ids, np.int64(lrg_id))
    subset = objects[left_ind:right_ind]
    subset[:,1:3] = np.radians(subset[:,1:3])
    total_neighbs+=len(subset)
    neighbs_unrot = np.ones((len(subset), 3))
    neighbs_unrot[:,1:3] = subset[:,1:3]
    mid = np.array([1., ra_mid, dec_mid], dtype = np.float64)
    right = np.array([1., ra_1, dec_1], dtype = np.float64)

    mid_rot, right_rot, neighbs_rot = rotate(mid, right, neighbs_unrot)
    
    lrg_radius = calc_distance((mid_rot[1], mid_rot[2]),
                               [(right_rot[1], right_rot[2])])[0] * lrg_add


    if (lrg_radius > 5 or lrg_radius < 3):
        oob_lrgs+=1
        continue
    #print lrg_radius

    r = calc_distance((mid_rot[1], mid_rot[2]), (neighbs_rot[:,1:]))
    r = np.multiply(r, lrg_add)
    oob_inds = np.where((r < -15) | (r > 15))[0]
    outOfGrid += len(oob_inds)
    r = np.delete(r, oob_inds)
    scaled_radii = np.divide(r, lrg_radius)
    x = np.multiply(map(np.cos, neighbs_rot[:,1]), scaled_radii)
    y = np.multiply(map(np.sin, neighbs_rot[:,1]), scaled_radii)
    
    xs = np.hstack((xs, x))
    ys = np.hstack((ys, y))
    
    xbs = map(lambda j: bs.bisect(xbin_ledge, j), x)
    ybs = map(lambda j: bs.bisect(ybin_ledge, j), y)
    def add_to_bins(ybin, xbin):
        grid[ybin, xbin] +=1
        return
    map(add_to_bins, ybs, xbs)
    
print '%d LRGs thrown out.' % oob_lrgs
print '%d total LRGs.' % len(pairs)
print '%d neighbors thrown out.' % outOfGrid
print '%d total neighbors.' % total_neighbs
#print grid

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
circ_1 = plt.Circle((1, 0), radius = .05, color = 'b')
circ_2 = plt.Circle((-1, 0), radius = .05, color = 'g')
circ_mid = plt.Circle((0, 0), radius = .05, color = 'r')

#plt.scatter(xs, ys)

x = np.linspace(-5, 5, xbins)
y = np.linspace(-5, 5, ybins)
X, Y = np.meshgrid(x, y)
plt.contour(x, y, grid)
#plt.scatter(xs, ys)
plt.xlabel('[R_sep]')
plt.ylabel('[R_sep]')
plt.title('Objects neighboring LRG pairs %d through %d'
          % (min(pair_range), max(pair_range)))
ax.add_artist(circ_1)
ax.add_artist(circ_2)
ax.add_artist(circ_mid)
#plt.savefig(out_dir + 'scatter_scaled_stacked_%d_%d'
#            % (min(pair_range), max(pair_range)))
plt.show()
