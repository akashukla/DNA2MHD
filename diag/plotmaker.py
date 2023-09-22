import dna2mhd_utils as dn
import output

lpaths = ["/scratch/08929/echansen/dna2mhdrun437"]
#fname = {lpaths[0]:"../input_files/DNAHD.out11263842",lpaths[1]:"../input_files/DNAHD.out11263842"}

opts = ['b','v','bdv','vdb','cbdb','bdcb','vdv','bdb']
for lpath in lpaths:
#    inds = output.maxinds(fname[lpath])
    inds = [[['V',26,10,63],'V'],[['V',3,27,50],'V'],[['V',7,36,35],'V'],[['V',1,37,38],'V'],[['V',10,19,1],'V'],[['V',30,38,2],'V']]
 #   dn.getb(lpath)
 #   dn.getv(lpath)
    print('Got NL Simulation')
    dn.plot_energy(lpath,2)
    dn.plot_enspec(lpath,4,log=True,newload=True)

    for ind in inds:
        for i in range(3):
            dn.plot_bv(lpath,ind[0][1],ind[0][2],ind[0][3],i,show=False,ask=False)
            dn.plot_nls(lpath,ind[0][1],ind[0][2],ind[0][3],i,show=False,ask=False)
            for opt in opts:
                dn.plot_bspectrum(lpath,ind[0][1],ind[0][2],ind[0][3],i,show=False,opt=opt)
    
        print("Plotted ",ind)
