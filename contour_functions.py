import numpy as np
import healpy as hp
import bisect as bs
import numpy.lib.recfunctions as rf
import matplotlib.pyplot as plt
import time

###LOCAL IMPORTS###
import cosmology_functions as cf
import spherical_geo as sg
from user_settings import user_settings as us
from contours_from_functions import contour_settings as cs

pair_adds = cs['USE_PAIR_ADDS']
scaled = cs['SCALED']
add_word = 'pair' if pair_adds else 'neighbor'
nside = cs['nside']

half_width = cs['max_distance'] if not scaled else cs['max_distance'] / 5.
half_height = cs['max_distance'] if not scaled else cs['max_distance'] / 5.

#Converts dec from the usual astrophysical convention to
#the otherwise usual convention, with theta=0 at the pole
#and theta varying from 0 to pi
def healpyDecConvention(decs):
    return np.multiply(np.subtract(decs, np.pi / 2.), -1.)

#Reads pairs and neighbors from the paths specified in contour_settings.py
#and returns structured arrays of the same catalogs, with ra/dec in
#radians and with pixel numbers added under the field 'pix'.
def addPixAndRadians():
    print 'Reading pairs...'
    pairs = np.genfromtxt(cs['pair_path'], delimiter=',', names=True)
    print 'Reading neighbors...'
    neighbors = np.genfromtxt(cs['neighbor_path'], delimiter=',', names=True)

    print 'Converting to radians...'
    for i in ['ra_mid', 'dec_mid', 'ra1', 'dec1', 'ra2', 'dec2']:
        pairs[i] = map(np.radians, pairs[i])
    for i in ['ra_gal', 'dec_gal']:
        neighbors[i] = map(np.radians, neighbors[i])

    print 'Calculating pixels...'
    pair_pix = hp.ang2pix(
        nside,
        healpyDecConvention(pairs['dec_mid']),
        pairs['ra_mid'],
        nest=True)
    neighbor_pix = hp.ang2pix(
        nside,
        healpyDecConvention(neighbors['dec_gal']),
        neighbors['ra_gal'],
        nest=True)
    pairs = rf.append_fields(pairs, 'pix', data=pair_pix, dtypes=int)
    neighbors = rf.append_fields(neighbors, 'pix', data=neighbor_pix, dtypes=int)
    return pairs, neighbors

def pixDict(neighbors):
    print 'Generating pixel dictionary...'
    pix_dict = {}
    pix_nums = neighbors['pix']
    #Add the index of each object as a value, with its pixel
    #number as key
    for n in range(len(neighbors)):
        pix_dict.setdefault(pix_nums[n], []).append(n)
    return pix_dict

def genContours(pairs, neighbors, pix_dict, chunknum, chunks):
    lb = (chunknum * len(pairs) / chunks)
    ub = ((chunknum+1) * len(pairs) / chunks)
    pairs = pairs[lb:ub]
    no_neighbors = 0
    oob_pairs = 0
    oob_neighbors = 0

    x_bins = np.linspace(-1. * half_width, half_width, cs['xbins'])[1:-1]
    y_bins = np.linspace(-1. * half_height, half_height, cs['ybins'])[1:-1]

    grid = np.zeros(shape=(cs['xbins'], cs['ybins']))

    t = time.time()
    print 'Adding to contours...'
    for p in range(len(pairs)):
        print (p + lb)
        curr = pairs[p]
        near_pix = np.append(hp.get_all_neighbours(nside, curr['pix']), curr['pix'])
        near_pix = near_pix[np.where(near_pix != -1)[0]]
        candidate_inds = []
        for j in near_pix:
            try:
                candidate_inds.append(pix_dict[j])
            except KeyError:
                pass
        candidate_inds = np.concatenate(candidate_inds)
        candidates = neighbors[np.array(candidate_inds)]
        pair_radius = cf.physicalDistance(
            (curr['ra_mid'], curr['dec_mid']),
            [(curr['ra1'], curr['dec1'])],
            [curr['z']]
            )
        distances = cf.physicalDistance(
            (curr['ra_mid'], curr['dec_mid']),
            zip(candidates['ra_gal'], candidates['dec_gal']),
            candidates['z'] if pair_adds else np.array([curr['z']] * len(candidates))
            )
        inbounds = candidates[np.where(distances < cs['max_distance'])[0]]
        oob_neighbors += len(candidates) - len(inbounds)
        distances = distances[np.where(distances < cs['max_distance'])[0]]
        unrotated_neighbors = np.ones((len(inbounds), 3))
        unrotated_neighbors[:,1] = inbounds['ra_gal']
        unrotated_neighbors[:,2] = inbounds['dec_gal']
        unrotated_mid = (1., curr['ra_mid'], curr['dec_mid'])
        unrotated_right = (1., curr['ra1'], curr['dec1'])
        rotated_mid, rotated_right, rotated_neighbors = sg.rotate(
            unrotated_mid,
            unrotated_right,
            unrotated_neighbors
            )
        if (pair_radius > 5 or pair_radius < 3):
            oob_pairs +=1
            continue
        if scaled:
            scaled_distances = np.divide(distances, pair_radius)
            xs = np.multiply(map(np.cos, rotated_neighbors[:,1]), scaled_distances)
            ys = np.multiply(map(np.sin, rotated_neighbors[:,1]), scaled_distances)
        else:
            xs = np.multiply(map(np.cos, rotated_neighbors[:,1]), distances)
            ys = np.multiply(map(np.sin, rotated_neighbors[:,1]), distances)
        xbs = map(lambda j: bs.bisect(x_bins, j), xs)
        ybs = map(lambda j: bs.bisect(y_bins, j), ys)
        def add_to_bins(ybin, xbin):
            grid[ybin, xbin] += 1
            return
        map(add_to_bins, ybs, xbs)

    return grid

def plotContours(grid):

    x = np.linspace(-1. * half_width, half_width, cs['xbins'])
    y = np.linspace(-1. * half_height, half_height, cs['ybins'])

    #print 'Saving grid...'
    #np.savetxt(cs['out_dir'] + 'grid_%s_adds_nside%d_%dx_%dy.csv' % (add_word, nside, cs['xbins'], cs['ybins']), grid, delimiter=',')
    plt.contour(x, y, grid)
    plt.savefig(cs['out_dir'] + 'plot_%s_adds_nside%d_%dx_%dy.png' % (add_word, nside, cs['xbins'], cs['ybins']))
    plt.show()
