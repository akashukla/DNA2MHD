import dna2mhd_utils as dn
import output

lpaths = ["/scratch/08929/echansen/dna2mhdrun374","/scratch/08929/echansen/dna2mhdrun375"]
fname = {lpaths[0]:"../input_files/DNAHD.out11263842",lpaths[1]:"../input_files/DNAHD.out11263842"}

for lpath in lpaths:
    inds = output.maxinds(fname[lpath])
    dn.getb(lpath)
    dn.getv(lpath)
    print('Got NL Simulation')
    dn.plot_energy(lpath,2)
    dn.plot_enspec(lpath,4,log=True,newload=True)

    for ind in inds:
        for i in range(3):
            dn.plot_bv(lpath,ind[0][1],ind[0][2],ind[0][3],i,show=False,ask=False)
            dn.plot_bspectrum(lpath,ind[0][1],ind[0][2],ind[0][3],i,show=False)
            dn.plot_vspectrum(lpath,ind[0][1],ind[0][2],ind[0][3],i,show=False)
        print("Plotted ",ind)