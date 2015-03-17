import numpy as np
import healpy as hp
import bisect as bs
import spherical_geo as sg
import cosmology_functions as cf
from healpy_settings import healpy_params as hpparam
import time
import sys
import os

n_path = '/home/chadavis/catalog_creation/LRG-pair-filaments/neighbors_with_pix.csv'
if hpparam['REAL_PAIRS']:
    p_path = '/home/chadavis/catalog_creation/LRG-pair-filaments/pairs_with_pix.csv'
else:
    p_path = '/home/chadavis/catalog_creation/LRG-pair-filaments/fake_cats/sampled_cat.csv'
out_dir = '/home/chadavis/catalog_creation/LRG-pair-filaments/matches/'

NSIDE = hpparam['NSIDE']

d = int(time.strftime("%d"))
m = int(time.strftime("%m"))
y = int(time.strftime("%y"))

chunks = 600
chunknum = int(os.environ['SGE_TASK_ID']) - 1
#chunknum = 5

print ''
print 'Chunk #%d' % chunknum

out_name = 'matches_%do%d_nside%d_%s_adds_%s.csv' % (chunknum, chunks, NSIDE, 'neighbor' if hpparam['USE_NEIGHBOR_ADDS'] else 'pair', 'real' if hpparam['REAL_PAIRS'] else 'fake')

print 'Reading neighbors...'
neighbors = np.genfromtxt(n_path, delimiter=',', dtype=str)
print 'Reading pairs...'
pairs = np.genfromtxt(p_path, delimiter=',', dtype=str)
gals = int(len(pairs) / chunks)
pairs = pairs[(chunknum * gals):((chunknum+1) * gals)]

print 'Pairs %d through %d' % ((chunknum * gals), (chunknum+1) * gals)

n_header = neighbors[0]
neighbors = neighbors[1:,]
p_header = pairs[0]
pairs = pairs[1:,]

nobjid = 0
pobjid = 7
lrgid = 0
ra_mid = 1
dec_mid = 2
pair_z = 3
ra_gal = 11
dec_gal = 12
z = 13


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

def distance(ra_p, dec_p, ra_n, dec_n, z_n):
    angle = sg.calc_distance(
        (np.radians(ra_p), np.radians(dec_p)),
        zip(
        map(np.radians, ra_n),
        map(np.radians, dec_n)
        )
        )
    adds = map(cf.angDiamDistSingle, z_n)
    return np.multiply(angle, adds)


pix_nums = neighbors[:,-1].astype(int)
print 'Creating pixel dictionary...'
pixDict = {}
for j in range(len(pix_nums)):
    pixDict.setdefault(pix_nums[j],[]).append(j)
    
matches = []

noneighbs = 0

print 'Matching...'
for p in pairs:
    pnum = int(p[-1])
    near_pix = hp.get_all_neighbours(NSIDE, pnum)
    inds = []
    for j in near_pix:
        try:
            inds.append(pixDict[int(j)])
        except KeyError:
            pass
    inds = np.concatenate(inds).astype(int)
    candidates = neighbors[inds]
    if hpparam['USE_NEIGHBOR_ADDS']:
        distances = distance(p[ra_mid].astype(float), p[dec_mid].astype(float),
                             candidates[:,ra_gal].astype(float),
                             candidates[:,dec_gal].astype(float),
                             candidates[:,z].astype(float))
    else:
        distances = distance(p[ra_mid].astype(float), p[dec_mid].astype(float),
                             candidates[:,ra_gal].astype(float),
                             candidates[:,dec_gal].astype(float),
                             np.array([float(p[pair_z])] * len(candidates)))
    inbounds = neighbors[inds[distances <= 15]]
    print '%d / %d' % (len(inbounds), len(inds))
    if len(inbounds) == 0:
        noneighbs+=1
        continue
    inbounds = inbounds[np.where(inbounds[:,nobjid] != p[pobjid])[0]]
    idnum_arr = np.array([p[lrgid]] * np.shape(inbounds)[0]).reshape(np.shape(inbounds)[0], 1)
    inbounds = np.hstack((inbounds, idnum_arr))
    matches.append(inbounds)

matches = np.concatenate(matches)
print 'Writing...'
np.savetxt(out_dir + out_name, matches, fmt='%s', delimiter=',')
print 'Written.'
print '%d w/o neighbors.' % noneighbs
