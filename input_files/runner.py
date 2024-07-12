from subprocess import call
import os

par = open(r'parameters','r')
a = par.readlines()
par.close()

TT = 5

Ns = [8,16,32]
dts = [0.01,0.02,0.05,0.1,0.2,0.5]
dtns = [1,2,3,4,5,6]
ios = [2,3,4,5]
r = 3

lstarts = ["nkx0 = ",
           "kxinit_max = ",
           "nky0 = ",
           "nkz0 = ",
           "kyinit_max = ",
           "kzinit_max = ",
           "diagdir = ",
           "intorder = ",
           "max_itime = ",
           "dt_max = "]

lineys = []
for j in lstarts:
    for ii,line in enumerate(a):
        if j in line:
            lineys.append(ii)
            continue
print(lineys)

runn = 0
for N in Ns:
    for dti in range(len(dtns)):
        for io in range(4):
            lpath = '/scratch/08929/echansen/dna2mhdrunctpush'+str(N)+str(ios[io])+str(dtns[dti])+str(r)
            for jj in range(len(lstarts)):
                if jj < 2:
                    a[lineys[jj]] = lstarts[jj] + str(N) + '\n'
                if jj >= 2 and jj < 6:
                    a[lineys[jj]] = lstarts[jj] + str(2*N) + '\n'
                if jj == 6:
                    a[lineys[jj]] = lstarts[jj] + "'" + lpath + "'" + '\n'
                if jj == 7:
                    a[lineys[jj]] = lstarts[jj] + str(ios[io]) + "\n"
                if jj == 8:
                    a[lineys[jj]] = lstarts[jj] + str(int(TT/dts[dti])+1) + '\n'
                if jj == 9:
                    a[lineys[jj]] = lstarts[jj] + str(dts[dti]) + '\n'
            par = open(r'parameters','w')
            par.writelines(a)
            par.close()
                        
            if os.path.exists(lpath):
                call(["rm","-rf",lpath])
            os.mkdir(lpath)
            tt = call(["ibrun","../bin2/dna"])
                        
            print("Done "+str(N)+str(ios[io])+str(dtns[dti])+str(r))
            runn += 1
