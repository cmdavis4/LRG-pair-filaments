import numpy as np
import healpy as hp
import sys

in_path = '/data2/scratch/pairs_6_10.txt'
out_dir = '/home/chadavis/catalog_creation/LRG-pair-filaments/'
out_filename='pairs_with_pix.csv'
res=16
ra_ind = 1
dec_ind = 2
dlim = ' '

print 'Reading neighbors...'
n = np.genfromtxt(in_path, dtype=str, delimiter=dlim)

print 'Calculating pixels...'
ras = np.radians(n[1:,ra_ind].astype(np.float64))
decs = (np.radians(n[1:,dec_ind].astype(np.float64)) - (np.pi / 2.)) * -1.
pix = hp.ang2pix(res, decs, ras)
pix = pix.reshape(len(pix), 1)
out = np.hstack((n[1:,], pix.astype(str)))
print 'Sorting...'
out = out[out[:,-1].argsort()]
header = np.append(n[0], 'pix')
out = np.vstack((header, out))

print 'Writing...'
np.savetxt(out_dir + out_filename, out, fmt='%s', delimiter=',')
