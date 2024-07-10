import os

os.system('module load remora')

par = open(r'parameters','r')
a = par.readlines()
par.close()

TT = 5

Ns = [8,16,32]
dts = [0.05,0.1,0.2]
dtns = [1,2,3]
ios = [2,3,4,5,41]
rs = [7,3]

lineys = [15,16,17,38,39,40,59,91,117,118,119,120,126,133,135]
lstarts = ["nkx0 = ",
           "nky0 = ",
           "nkz0 = ",
           "kxinit_max = ",
           "kyinit_max = ",
           "kzinit_max = ",
           "diagdir = ",
           "intorder = ",
           "en_leftwhist = ",
           "en_leftcyclo = ",
           "en_rightwhist = ",
           "en_rightcyclo = ",
           "random_state = ",
           "max_itime = ",
           "dt_max = "]

for N in Ns:
    for dti in range(dtns):
        for i in range(5):
            for r in rs:
                lpath = 'dna2mhdrunct'+str(N)+str(ios[i])+str(dtns[dti])+str(r)
                for i in range(len(lineys)):
                    if (i == 0) or (i == 3):
                        a[lineys[i]] = lstarts[i] + str(N) + '\n'
                    elif (i < 6):
                        a[lineys[i]] = lstarts[i] + str(2*N) + '\n'
                    elif i == 6:
                        a[lineys[i]] = lstarts[i] + lpath + '\n'
                    elif i == 7:
                        a[lineys[i]] = lstarts[i] + ios[i] + '\n'
                    elif i < 12:
                        if r == 7:
                            a[lineys[i]] = lstarts[i] + str(0.25) + '\n'
                        elif i == 9:
                            a[lineys[i]] = lstarts[i] + str(1.0) + '\n'
                        else:
                            a[lineys[i]] = lstarts[i] + str(0.0) + '\n'
                    elif i == 12:
                        a[lineys[i]] = lstarts[i] + str(r) + '\n'
                    elif i == 13:
                        if ios[i] == "41":
                            a[lineys[i]] = lstarts[i] + str(int(50/dts[dti])+1) + '\n'
                        else:
                            a[lineys[i]] = lstarts[i] + str(int(TT/dts[dti])+1) + '\n'
                    else:
                        a[lineys[i]] = lstarts[i] + str(dts[dti]) + '\n'
                par = open(r'parameters','w')
                for i in a:
                    par.write(i)
                par.close()
                os.system('mkdir '+lpath)
                os.system('remora ibrun ../bin2/dna')
                quit
                    
                
                

