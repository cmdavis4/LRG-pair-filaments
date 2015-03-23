import contour_functions as cf
from numpy import savez, load
from contours_from_functions import contour_settings as cs

grid = load('/home/chadavis/catalog_creation/LRG-pair-filaments/steps/temp_files/step2_0_o_%d_%s_%s_adds_%dx%d_nside%d.npz' % (cs['chunks'], real, add, xb, yb, nside), grid=grid)['grid']

for j in range(0, cs['chunks'] + 1):
    grid += load('/home/chadavis/catalog_creation/LRG-pair-filaments/steps/temp_files/step2_%d_o_%d_%s_%s_adds_%dx%d_nside%d.npz' % (j, cs['chunks'], real, add, xb, yb, nside), grid=grid)['grid']

savez('/home/chadavis/catalog_creation/LRG-pair-filaments/steps/temp_files/step2_complete_%s_%s_adds_%dx%d_nside%d.npz' % (real, add, xb, yb, nside), grid)

cf.plotContours(grid)
