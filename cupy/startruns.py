from subprocess import run

wave = False
nlinv = True
re = True
if wave:
    ii = [0,6,11,17]

    if not re:
        for i in ii:
            run(["sbatch","-J","wave"+str(i),"submitwave.cmd",str(i)])
    else:
        for i in ii:
            run(["sbatch","-J","2wave"+str(i),"resubmitwave.cmd",str(i)])

if nlinv:
    ii = range(0,6)

    if not re:
        for i in ii:
            run(["sbatch","-J","nlinv"+str(i),"submitnlinv.cmd",str(i)])
    else:
        for i in ii:
            run(["sbatch","-J","2nlinv"+str(i),"resubmitnlinv.cmd",str(i)])
