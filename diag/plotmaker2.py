import dna2mhd_utils as dn
import output

lpaths = ["/scratch/08929/echansen/dna2mhdrun378"]
fname = {lpaths[0]:"../input_files/DNAHD.out11320789"}
inds = output.maxinds(fname["/scratch/08929/echansen/dna2mhdrun378"])
for lpath in lpaths:
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
