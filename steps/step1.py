import contour_functions as cf
from numpy import savez

pairs, neighbors = cf.addPixAndRadians()
pix_dict = cf.pixDict(neighbors)
savez('/home/chadavis/catalog_creation/LRG-pair-filaments/steps/temp_files/step1.npz', pairs=pairs, neighbors=neighbors, pix_dict=pix_dict)
