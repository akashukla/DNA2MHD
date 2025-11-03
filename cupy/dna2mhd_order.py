from dna2mhd_cupy import DNA2MHD

N = 288
nkx0 = N
nky0 = N
nkz0 = 128

kxmin = 0.025
kymin = 0.025
kzmin = 0.005

nu = 0.0
eta = 0.0
dt = 0.1
iterations = 10
lpath = "/pscratch/sd/e/echansen/testingcupy1"

solver = DNA2MHD(nkx0,nky0,nkz0,kxmin,kymin,kzmin,
                 nu,eta,
                 dt,iterations,
                 lpath,linear=False,explicitrk4=True,
                 initialcondition="hallwave",energystart=0.01,init_kolm=0,hmhdwave=[1,0,0,0],
                 forcetype="hallwave",forceamp=0.0,nforce=4,forcewave=[1,0,0,0],
                 hyper=1,hallparam=1.0,
                 solveprec=16,maxwallclock=86200)
solver.splitsimulation()
