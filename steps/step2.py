import contour_functions as cf
from numpy import savez, load
import os
from contours_from_functions import contour_settings as cs

inp = load('/home/chadavis/catalog_creation/LRG-pair-filaments/steps/temp_files/step1.npz')

pairs = inp['pairs']
neighbors = inp['neighbors']
pix_dict=inp['pix_dict']

chunknum = int(os.environ['SGE_TASK_ID']) - 1

grid = cf.genContours(pairs, neighbors, pix_dict, chunknum, cs['chunks'])

add = 'pair' if cf['USE_PAIR_ADDS'] else 'neighbor'
real = 'real' if cf['REAL_PAIRS'] else 'fake'
xb = cf['xbins']
yb = cf['ybins']
nside = cf['nside']

savez('/home/chadavis/catalog_creation/LRG-pair-filaments/steps/temp_files/step2_%d_o_%d_%s_%s_adds_%dx%d_nside%d.npz' % (chunknum, cs['chunks'], real, add, xb, yb, nside), grid=grid)
