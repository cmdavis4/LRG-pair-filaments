#+++++++++++++++++++++++++++++++++++++++
#
# Cosmology functions
#
# Description:
# This file contains the functions
# for all cosmology calculations
# 
#---------------------------------------

# local imports
from cosmology_settings import cosmology_params as cosmo_params

def angDiamDistSingle(z):
    """Angular diameter distance of an object at redshift z"""
    E_inv = lambda z: 1. / math.sqrt((cosmo_param['omega_m'] * (1. + z)**3) + (cosmo_param['omega_k'] * (1+z)**2) + cosmo_param['omega_l'])
    Dc = cosmo_param['Dh'] * si.quad(E_inv, 0, z)[0]
    Dm = Dc
    Da = (Dm / (1. + z)) * (1. / .7)
    return Da


def angDiamDist(z1, z2):
    """Angular diameter distance of an object at redshift z2 w.r.t. an
observer at redshift z1"""
    E_inv = lambda z: 1 / math.sqrt((cosmo_param['omega_m'] * (1 + z)**3) + (cosmo_param['omega_k'] * (1+z)**2) + cosmo_param['omega_l'])
    Dm = lambda z: cosmo_param['Dh'] * si.quad(E_inv, 0, z)[0]
    Da12 = (1 / (1 + z2)) * (Dm(z2) - Dm(z1))
    return Da12
