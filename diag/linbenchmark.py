import dna2mhd_utils as dn
import numpy as np

lpath = '/scratch/08929/echansen/dna2mhdrun308'

inds = []
for i in range(64):
    for j in range(64):
        for k in range(32):
            inds.append([i,j,k])

chosen = []
for l in range(50):
    chosen.append(inds.pop(np.random.random_integers(0,64*64*32-1)))

for arr in chosen:
    dn.plot_vspectrum(lpath,arr[0],arr[1],arr[2],0)
