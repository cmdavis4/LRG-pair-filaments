import numpy as np
import healpy as hp
import bisect as bs
import spherical_geo as sg
import cosmology_functions as cf
from healpy_settings import healpy_params as hpparam

n_path = '/home/chadavis/catalog_creation/LRG-pair-filaments/neighbors_with_pix.csv'
p_path = '/home/chadavis/catalog_creation/LRG-pair-filaments/pairs_with_pix.csv'
out_dir = '/home/chadavis/catalog_creation/LRG-pair-filaments/'
out_name = 'matches.csv'

print 'Reading neighbors...'
neighbors = np.genfromtxt(n_path, delimiter=',', dtype=str)
print 'Reading pairs...'
pairs = np.genfromtxt(p_path, delimiter=',', dtype=str)

n_header = neighbors[0]
neighbors = neighbors[1:,]
p_header = pairs[0]
pairs = pairs[1:,]

nobjid = 0
pobjid = 7
lrgid = 0
ra_mid = 1
dec_mid = 2
ra_gal = 11
dec_gal = 12
z = 13

NSIDE = hpparam['NSIDE']

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
    lind = lambda x: bs.bisect_left(neighbors[:,-1].astype(int), x)
    rind = lambda x: bs.bisect_right(neighbors[:,-1].astype(int), x)
    get_inds = lambda x: range(lind(x), rind(x))
    inds = np.concatenate(map(get_inds, p))
    return np.unique(inds).astype(int)

def distance(ra_p, dec_p, ra_n, dec_n, z_n):
    angle = sg.calc_distance((ra_p, dec_p), zip(ra_n, dec_n))
    adds = map(cf.angDiamDistSingle, z_n)
    return np.multiply(angle, adds)


pix_nums = neighbors[:,-1].astype(int)
print 'Creating pixel dictionary...'
pixDict = {}
for j in range(len(pix_nums)):
    pixDict.setdefault(pix_nums[j],[]).append(j)
    
matches = []

print 'Matching...'
for p in pairs:
    pnum = int(p[-1])
    near_pix = hp.get_all_neighbours(NSIDE, pnum)
    print pnum
    inds = []
    for j in near_pix:
        try:
            inds.append(pixDict[int(j)])
        except KeyError:
            pass
    inds = np.concatenate(inds).astype(int)
    candidates = neighbors[inds]
    print np.shape(candidates)
    distances = distance(p[ra_mid].astype(float), p[dec_mid].astype(float),
                         candidates[:,ra_gal].astype(float),
                         candidates[:,dec_gal].astype(float),
                         candidates[:,z].astype(float))
    inbounds = neighbors[inds[distances <= 15]]
    inbounds = inbounds[np.where(inbounds[:,nobjid] != p[pobjid])[0]]
    idnum_arr = np.array([p[lrgid]] * np.shape(inbounds)[0]).reshape(np.shape(inbounds)[0], 1)
    inbounds = np.hstack((inbounds, idnum_arr))
    matches.append(inbounds)
    
    
np.savetxt(out_dir + out_name, matches, fmt='%s', delimiter=',')
