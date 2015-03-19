contour_settings = {
    'USE_PAIR_ADDS':True,
    'REAL_PAIRS':True,
    'SCALED':False,
    'xbins':20,
    'ybins':20,
    'out_name':'test_fig.png',
    'max_distance':15.,
    'nside':16,
    'out_dir':'/home/chadavis/catalog_creation/LRG-pair-filaments/results/',
    'pair_path':'/home/chadavis/catalog_creation/LRG-pair-filaments/pairs_6_10.csv',
    'neighbor_path':'/home/chadavis/catalog_creation/LRG-pair-filaments/dr10_allspec.csv'
    }

if __name__ == "__main__":
    import contour_functions as cf
    pairs, neighbors = cf.addPixAndRadians()
    pd = cf.pixDict(neighbors)
    grid = cf.genContours(pairs, neighbors, pd, 5, 600)
    cf.plotContours(grid)

