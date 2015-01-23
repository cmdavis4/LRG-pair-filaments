import numpy as np
import matplotlib.pyplot as plt
import bisect as bs
import math
import scipy.integrate as si
import numpy.linalg as LA

pair_range = range(10000, 11000)
xbins = 10
ybins = 10
xinc = 6. / xbins
yinc = 6. / ybins
xbin_ledge = np.linspace(-3, 3, xbins)
ybin_ledge = np.linspace(-3, 3, ybins)
neighbor_id_ind = 0
neighbor_ra_ind = 1
neighbor_dec_ind = 2
neighbor_z_ind = 3
lrg_id_ind = 0
omega_m = .3
omega_l = .7
omega_k = 0
Dh = 3000. #h^-1 * Mpc
neighbor_path = 'neighbors.cat'
pair_path = 'pairs.cat'
out_dir = ''
#neighbor_path = '/home/chadavis/catalog_creation/astro_image_processing/LRG/data/ra_dec_z.cat'
#pair_path = '/data2/scratch/pairs_6_10.txt'
neighbors = []

#Angular diameter distance of an object at redshift z
def angDiamDistSingle(z):
    E_inv = lambda z: 1. / math.sqrt((omega_m * (1. + z)**3) + (omega_k * (1+z)**2) + omega_l)
    Dc = Dh * si.quad(E_inv, 0, z)[0]
    Dm = Dc
    Da = (Dm / (1. + z)) * (1. / .7)
    return Da

#Angular diameter distance of an object at redshift z2 w.r.t. an
#observer at redshift z1
def angDiamDist(z1, z2):
    E_inv = lambda z: 1 / math.sqrt((omega_m * (1 + z)**3) + (omega_k * (1+z)**2) + omega_l)
    Dm = lambda z: Dh * si.quad(E_inv, 0, z)[0]
    Da12 = (1 / (1 + z2)) * (Dm(z2) - Dm(z1))
    return Da12

def to_sphere(point):
    """accepts a x,y,z, coordinate and returns the same point in spherical coords (r, theta, phi) i.e. (r, azimuthal angle, altitude angle) with phi=0 in the zhat direction and theta=0 in the xhat direction"""
    coord = (np.sqrt(np.sum(point**2)),
             np.arctan2(point[1],point[0]),
             np.pi/2.0 - np.arctan2(np.sqrt(point[0]**2 + point[1]**2),point[2])
             )
    
    return np.array(coord)
def to_cartesian(point):
    """accepts point in spherical coords (r, theta, phi) i.e. (r, azimuthal angle, altitude angle) with phi=0 in the zhat direction and theta=0 in the xhat direction and returns a x,y,z, coordinate"""
    coord = point[0]*np.array([np.cos(point[2])*np.cos(point[1]),
                                 np.cos(point[2])*np.sin(point[1]),
                                 np.sin(point[2])])
    
    
    return coord
def rotate(p1, p2, p3):
    """rotates coordinate axis so that p1 is at the pole of a new spherical coordinate system and p2 lies on phi (or azimuthal angle) = 0

inputs:
p1: vector in spherical coords (phi, theta, r) where phi is azimuthal angle (0 to 2 pi), theta is zenith angle or altitude (0 to pi), r is radius
p2: vector of same format 
p3: vector of same format

Output:
s1:vector in (r, theta, phi) with p1 on the z axis
s2:vector in (r,, theta, phi) with p2 on the phi-hat axis
s3:transformed vector of p3
"""

    p1_cc=to_cartesian(p1) 
    p2_cc=to_cartesian(p2) 
    p3_cc=to_cartesian(p3) 

    p1norm = p1_cc/LA.norm(p1_cc)
    p2norm = p2_cc/LA.norm(p2_cc)
    p3norm = p3_cc/LA.norm(p3_cc)
    
    zhat_new =  p1norm
    x_new = p2norm - np.dot(p2norm, p1norm) * p1norm
    xhat_new = x_new/LA.norm(x_new)
    
    yhat_new = np.cross(zhat_new, xhat_new)
    
    s1 = np.array(map(lambda x: np.dot(x, p1_cc), (xhat_new, yhat_new, zhat_new)))
    s2 = np.array(map(lambda x: np.dot(x, p2_cc), (xhat_new, yhat_new, zhat_new)))
    s3 = np.array(map(lambda x: np.dot(x, p3_cc), (xhat_new, yhat_new, zhat_new)))

    s1=to_sphere(s1) 
    s2=to_sphere(s2) 
    s3=to_sphere(s3)
    
    return s1, s2, s3


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

###############

print('Reading neighbors...')
neighbors = np.loadtxt(neighbor_path)
ids = np.array(neighbors[:,0], dtype = np.int64)        

print('Reading pairs...')
pairs = np.loadtxt(pair_path, skiprows = 1)

#####DELETE AFTER TESTING#########
pair_range = range(len(pairs))
##################################

objects = np.zeros(shape=(len(neighbors), 4))
objects[:,0] = neighbors[:,neighbor_id_ind] #id
objects[:,1] = neighbors[:,neighbor_ra_ind] #ra
objects[:,2] = neighbors[:,neighbor_dec_ind] #dec
objects[:,3] = neighbors[:,neighbor_z_ind] #z

print('Sorting neighbors...')
objects = np.array(sorted(objects, key = lambda x: x[0]), dtype=np.float64)
ids = sorted(ids)

offset = 0
grid = np.zeros(shape=(xbins, ybins))
outOfGrid = 0
total_neighbs = 0
oob_lrgs = 0
for i in pair_range:
    curr = pairs[i]
    lrg_id = curr[lrg_id_ind]

    ra_mid = np.radians(curr[1])
    dec_mid = np.radians(curr[2])
    ra_1 = np.radians(curr[5])
    dec_1 = np.radians(curr[6])
    ra_2 = np.radians(curr[10])
    dec_2 = np.radians(curr[11])
    lrg_add = angDiamDistSingle(curr[3])
    
    left_ind = bs.bisect_left(ids, np.int64(lrg_id))
    right_ind = bs.bisect(ids, np.int64(lrg_id))
    subset = objects[left_ind:right_ind]
    subset[:,1:3] = np.radians(subset[:,1:3])
    total_neighbs+=len(subset)
    adds = np.zeros(len(subset))
    neighbs_rot = np.zeros((len(subset), 3))
    mid = np.array([1., ra_mid, dec_mid], dtype = np.float64)
    right = np.array([1., ra_1, dec_1], dtype = np.float64)

    for j in range(len(subset)):
        neighb = np.array([1., subset[j, 1], subset[j, 2]], dtype = np.float64)
        mid_rot, right_rot, neighbs_rot[j,] = rotate(mid, right, neighb)
        adds[j] = angDiamDistSingle(subset[j,3])
    
        
    lrg_radius = calc_distance(mid_rot[1], mid_rot[2],
                               right_rot[1], right_rot[2]) * lrg_add

    
    print (calc_distance(mid[1], mid[2],
                        right[1], right[2]) * lrg_add)
    print lrg_radius
    

    if (lrg_radius > 10 or lrg_radius < 6):
        oob_lrgs+=1
        continue
    #print lrg_radius

    for j in range(len(subset)):
        r = calc_distance(mid_rot[1], mid_rot[2],
                          neighbs_rot[j, 1], neighbs_rot[j, 2]) * lrg_add
        if (r < -15 or r > 15):
            outOfGrid += 1
            continue
        #print r
        scaled_radius = r / lrg_radius
        x = np.cos(neighbs_rot[j, 1]) * scaled_radius
        y = np.sin(neighbs_rot[j, 1]) * scaled_radius
        #print '%f, %f' % (x, y)
        #print ''
        xbin = bs.bisect(xbin_ledge, x) - 1
        ybin = bs.bisect(ybin_ledge, y) - 1
        grid[xbin, ybin] += 1

print '%d LRGs thrown out.' % oob_lrgs
print '%d total LRGs.' % len(pairs)
print '%d neighbors thrown out.' % outOfGrid
print '%d total neighbors.' % total_neighbs

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
circ_1 = plt.Circle((.5, 0), radius = .05, color = 'b')
circ_2 = plt.Circle((-.5, 0), radius = .05, color = 'g')
circ_mid = plt.Circle((0, 0), radius = .05, color = 'r')
x = np.linspace(-3, 3, xbins)
y = np.linspace(-3, 3, ybins)
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
plt.savefig(out_dir + '/scatter_scaled_stacked_%d_%d'
            % (min(pair_range), max(pair_range)))
plt.show()
