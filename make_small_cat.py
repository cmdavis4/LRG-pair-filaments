import numpy as np
import random
import bisect as bs

neighbor_path = '/home/chadavis/catalog_creation/astro_image_processing/LRG/data/ra_dec_z.cat'
pair_path = '/data2/scratch/pairs_6_10.txt'
out_path = '/home/chadavis/catalog_creation/'
neighbor_ra_ind = 2
neighbor_dec_ind = 3
num_pairs = 10

neighbors = np.loadtxt(neighbor_path)
ids = np.array(neighbors[:,0], dtype = np.int64)        

pairs = np.loadtxt(pair_path, skiprows = 1)

objects = np.zeros(shape=(len(neighbors), 4))
objects[:,0] = neighbors[:,0] #id
objects[:,1] = neighbors[:,neighbor_ra_ind] #ra
objects[:,2] = neighbors[:,neighbor_dec_ind] #dec
objects[:,3] = neighbors[:,4] #z

objects = np.array(sorted(objects, key = lambda x: x[0]), dtype=np.float64)
ids = sorted(ids)

p_subset = []
n_subset = np.array([]).reshape(0,len(objects[0]))
for i in range(num_pairs):
    p_subset.append(pairs[random.randint(0, len(pairs) - 1)])
p_subset = np.array(p_subset)

for i in range(len(p_subset)):
    #print p_subset[i,0]
    l_ind = bs.bisect_left(ids, p_subset[i, 0])
    #print(ids[l_ind])
    #print(ids[l_ind + 1])
    r_ind = bs.bisect(ids, p_subset[i, 0])
    print p_subset[i, 0]
    n_subset = np.vstack((n_subset, objects[l_ind:r_ind,]))
    print(n_subset)
    print ''

print n_subset
np.savetxt(out_path + 'neighbors.cat', n_subset)
np.savetxt(out_path + 'pairs.cat', p_subset)
