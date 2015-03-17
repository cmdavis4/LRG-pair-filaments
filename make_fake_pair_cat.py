#! /usr/global/paper/bin/python
from math import pi
import numpy as np
from scipy import integrate
from scipy.interpolate import interp1d
import pyfits as pf
import sys
import os

ndiv = 70
div = int(os.environ['SGE_TASK_ID']) - 1
#div = 0

indir = '/home/clampitt/filaments/'
outdir = '/home/chadavis/catalog_creation/LRG-pair-filaments/fake_cats/'

dr = 100.

# Radial binning
Rmin = 0.1
Rmax = 24.                              # Mpc/h

outfile = 'pair_cat_LRG_Rmax%.1f_rlos%.1f_%do%d.dat' % (Rmax, dr, div, ndiv)

################
def comoving_dist(z):     # kpc/h
    om_m = 0.3
    om_L = 0.7
    c = 3e5               # km/s
    H0 = 100.             # km h / (s Mpc)
    def dist_int(z):
        return 1./np.sqrt(om_m*(1.+z)**3 + om_L)
    res, err = integrate.quad(dist_int, 0., z)
    return [res * (c/H0) * 1000., err * (c/H0) * 1000.]
################

# Read lens coordinates
lensfile = 'LRG_DR7-Full.fits'

hdu = pf.open(indir +lensfile)
header = hdu[1].header
columns = hdu[1].columns
keys = columns.names

print 'columns = ', columns, '\n'
#print 'keys = ', keys, '\n'
print 'num gals = ', len(hdu[1].data.field(keys[0]))

ra_in = hdu[1].data.field('RA')
dec_in = hdu[1].data.field('DEC')
dec = dec_in * pi / 180.
z = hdu[1].data.field('Z')

nlens = len(ra_in)
print nlens, ' lenses'

# Work in proper distance along los
z_int = np.linspace(0.001, 0.81, 2000)
chi_int = []
for j in range(len(z_int)):
    res, err = comoving_dist(z_int[j])
    chi_int.append(res)
chi_int = np.array(chi_int)                         # kpc/h

func = interp1d(z_int, chi_int, kind='linear')
chi = func(z)
# los distance and ang. diameter distance same for no curvature
rlos = chi / (1.+z) / 1000.                 # Mpc/h

# Decide which galaxies to do
gal_id = np.arange(nlens)

div_edge = np.linspace(gal_id[0], gal_id[-1], ndiv+1)
this_id = np.where((div_edge[div] < gal_id) & (gal_id < div_edge[div+1]))[0]

# Center on lens
trigger = 0
print 'iterating over lenses'
for i in this_id:

    print 'lens number \t ', i, '\t', i/float(nlens)*100, '%'

    # Select nearby lenses
    tmax = Rmax / rlos[i]

    ra = ra_in * np.cos(dec[i]) * pi / 180.
    dist = np.sqrt((ra - ra[i])**2 + (dec - dec[i])**2)
    
    rcon = (dist < tmax)

    # Choose distant lenses in LOS separation
    loscon = ( (((rlos[i] + dr) < rlos) & (rlos < (rlos[i] + dr + 20.))) |
               (((rlos[i] - dr - 20.) < rlos) & (rlos < (rlos[i] - dr))) )
    print np.sum(loscon)
    
    nbr = np.where(rcon & loscon)[0]
    if (len(nbr) < 2): continue

    itself = np.where(nbr == i)[0]
    nbr = np.delete(nbr, [itself])

    # Find physical distance
    x_near = (ra[nbr] - ra[i]) * (rlos[i]+rlos[nbr])/2.
    y_near = (dec[nbr] - dec[i]) * (rlos[i]+rlos[nbr])/2.

    R_near = np.sqrt(x_near**2 + y_near**2)

    drlos = (rlos[i] - rlos[nbr])
    dz = (z[i] - z[nbr])

    ct = 0
    for j in nbr:
        dec_mid = (dec_in[i] + dec_in[j]) / 2.
        ra_mid = (ra[i] + ra[j]) / 2. / np.cos(dec_in[i] * pi/180.) * 180./pi
        row = np.hstack((ra_mid, dec_mid, (z[i]+z[j])/2., i, ra_in[i], dec_in[i], z[i], 
                         j, ra_in[j], dec_in[j], z[j], R_near[ct], drlos[ct], dz[ct]))
        ct = ct + 1

        if (j < i): continue                 # Assure each pair only counted once
        
        if (trigger == 0): tab = row
        else: tab = np.vstack((tab, row))

        trigger = 1

np.savetxt(outdir +outfile, tab, '%.4f')
