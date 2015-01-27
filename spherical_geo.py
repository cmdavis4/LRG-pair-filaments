#+++++++++++++++++++++++++++++++++++++++
#
# Spherical Geometry
#
# Description:
# This file contains the functions for
# spherical geometry calcuations 
# including rotations, coordinate transformations
# and angular separations
# 
#---------------------------------------

import numpy as np

def to_sphere(point):
    """accepts a x,y,z, coordinate and returns the same point in spherical coords (r, theta, phi) i.e. (r, azimuthal angle, altitude angle) with phi=0 in the zhat direction and theta=0 in the xhat direction"""
    coord = (np.sqrt(np.sum(point**2)),
             np.arctan2(point[1],point[0]),
             np.pi/2.0 - np.arctan2(np.sqrt(point[0]**2 + point[1]**2),point[2])
             )
    
    return np.array(coord)
def to_cartesian(point):
    """accepts point in spherical coords (r, theta, phi) i.e. (r, azimuthal angle, altitude angle) with phi=0 in the zhat direction and theta=0 in the xhat direction and returns a x,y,z, coordinate"""
    coord = point[0]*np.array([np.cos(point[2])*np.cos(point[1]),
                                 np.cos(point[2])*np.sin(point[1]),
                                 np.sin(point[2])])
    
    
    return coord
def rotate(p1, p2, p3):
    """rotates coordinate axis so that p1 is at the pole of a new spherical coordinate system and p2 lies on phi (or azimuthal angle) = 0

inputs:
p1: vector in spherical coords (phi, theta, r) where phi is azimuthal angle (0 to 2 pi), theta is zenith angle or altitude (0 to pi), r is radius
p2: vector of same format 
p3: vector of same format

Output:
s1:vector in (r, theta, phi) with p1 on the z axis
s2:vector in (r,, theta, phi) with p2 on the phi-hat axis
s3:transformed vector of p3
"""

    p1_cc=to_cartesian(p1) 
    p2_cc=to_cartesian(p2) 
    p3_cc=to_cartesian(p3) 

    p1norm = p1_cc/LA.norm(p1_cc)
    p2norm = p2_cc/LA.norm(p2_cc)
    p3norm = p3_cc/LA.norm(p3_cc)
    
    zhat_new =  p1norm
    x_new = p2norm - np.dot(p2norm, p1norm) * p1norm
    xhat_new = x_new/LA.norm(x_new)
    
    yhat_new = np.cross(zhat_new, xhat_new)
    
    s1 = np.array(map(lambda x: np.dot(x, p1_cc), (xhat_new, yhat_new, zhat_new)))
    s2 = np.array(map(lambda x: np.dot(x, p2_cc), (xhat_new, yhat_new, zhat_new)))
    s3 = np.array(map(lambda x: np.dot(x, p3_cc), (xhat_new, yhat_new, zhat_new)))

    s1=to_sphere(s1) 
    s2=to_sphere(s2) 
    s3=to_sphere(s3)
    
    return s1, s2, s3


def calc_distance(ra1, dec1, ra2, dec2):
    '''Calculate the circular angular distance of two points on a sphere.'''
    lambda_diff = ra1  - ra2
    cos1 = np.cos(dec1)
    cos2 = np.cos(dec2)
    sin1 = np.sin(dec1)
    sin2 = np.sin(dec2)
    
    num = (cos2 * np.sin(lambda_diff)) ** 2.0 + (cos1 * sin2 - sin1 * cos2 * np.cos(lambda_diff)) ** 2.0
    denom = sin1 * sin2 + cos1 * cos2 * np.cos(lambda_diff)
    
    return np.arctan2(np.sqrt(num), denom)
