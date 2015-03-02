import numpy as np
import healpy as hp
import bisect as bs
import spherical_geo as sg
import cosmology_functions as cf

n_path = ''
p_path = ''
out_dir = ''
out_name = ''

neighbors = np.genfromtxt(n_path, delimiter=',', dtype=str, names=True)
pairs = np.genfromtxt(p_path, delimiter=',', dtype=str, names=True)
NSIDE = 16

'''
Steps for one pair:
-Read in its pixel number
-Get the numbers of the nearest pixels
-Find all objects in the neighbors catalog with those pixel numbers
-Calculate the distance from the midpoint of the LRG pair for each
 of these objects, using the ADD of the neighbor, not of the pair
-Objects within 15 Mpc are neighbors
'''

matches = np.array([])

def matchPix(p):
    lind = lambda x: bs.bisect_left(pairs['pix'].astype(int), x)
    rind = lambda x: bs.bisect_right(pairs['pix'].astype(int), x)
    return zip(map(lind, p), map(rind, p))

def distance(ra_p, dec_p, ra_n, dec_n, z_n):
    angle = sg.calc_distance(ra_p, dec_p, ra_n, dec_n)
    adds = map(cf.angDiamDistSingle, z_n)
    return np.multiply(angle, adds)

for n in neighbors:
    pnum = float(n['pix'])
    near_pix = hp.get_all_neighbours(NSIDE, pnum)
    inds = matchPix(p)
    candidates = neighbors[inds]
    distances = distance(n['ra_mid'].astype(float), n['dec_mid'].astype(float),
                         candidates['ra_gal'].astype(float),
                         candidates['dec_gal'].astype(float),
                         candidates['z'].astype(float))
    inbounds = neighbors[inds[distances <= 15]]
    inbounds = hstack((inbounds, [n['id']] * len(inbounds)))
    matches = vstack((matches, inbounds))
    
np.savetxt(out_dir + out_name, matches, fmt='%s', delimiter=',')
    