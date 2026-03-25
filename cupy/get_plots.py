from dna2mhd_exp import DNA2MHD
import dna2mhd_utils as dn
import numpy as np
import os


# Specific to 5 interaction, 5 setting file!
for seti in range(0,5):
    for intnum in range(0,5):

        with h5py.File("initconds032426.hdf5","w") as f:
            fname = f["fnames"][seti]

            lpath = "/pscratch/sd/e/echansen/threewaves032426/"+fname+str(intnum)+"/"
        dn.plot_energy(lpath)
        dn.plot_enspec(lpath,zz=-1,version=0)
        dn.mode_break(lpath,show=False)
        dn.threewaveenergy(lpath)
        dn.nonlinearities(lpath)
