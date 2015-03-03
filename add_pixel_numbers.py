import numpy as np
import healpy as hp
import sys

#in_path = '/data3/LRGPairs/dr10_allspec.csv'
#out_filename='neighbors_with_pix.csv'
in_path = '/data2/scratch/pairs_6_10.txt'
out_filename='pairs_with_pix.csv'
out_dir = '/home/chadavis/catalog_creation/LRG-pair-filaments/'
res=16

#####
#ra_ind = 11
#dec_ind = 12
#dlim = ','
ra_ind = 1
dec_ind = 2
#####

print 'Reading neighbors...'
n = np.genfromtxt(in_path, dtype=str)
header = n[0]
header = np.append(header, 'pix')
n = n[1:]
print n[0]

print 'Calculating pixels...'
ras = np.radians(n[:,ra_ind].astype(np.float64))
decs = (np.radians(n[:,dec_ind].astype(np.float64)) - (np.pi / 2.)) * -1.
pix = hp.ang2pix(res, decs, ras)
pix = pix.reshape(len(pix), 1)
out = np.hstack((n, pix.astype(str)))
print 'Sorting...'
out = out[out[:,-1].astype(int).argsort()]
out = np.vstack((header, out))

print 'Writing...'
np.savetxt(out_dir + out_filename, out, fmt='%s', delimiter=',')
