import numpy as np
import healpy as hp
import sys

neighbors_path = '/data3/LRGPairs/dr10_allspec.csv'
out_dir = '/home/chadavis/catalog_creation/LRG-pair-filaments/'
out_filename='neighbors_with_pix.csv'
res=16

print 'Reading neighbors...'
n = np.genfromtxt(neighbors_path, delimiter=',', dtype=str)

print 'Calculating pixels...'
ras = np.radians(n[1:,11].astype(np.float64))
decs = (np.radians(n[1:,12].astype(np.float64)) - (np.pi / 2.)) * -1.
pix = hp.ang2pix(res, decs, ras)
pix = pix.reshape(len(pix), 1)
out = np.hstack((n[1:,], pix.astype(str)))
print 'Sorting...'
out = out[out[:,-1].argsort()]
header = np.append(n[0], 'pix')
out = np.vstack((header, out))

print 'Writing...'
np.savetxt(out_dir + out_filename, out, fmt='%s', delimiter=',')
