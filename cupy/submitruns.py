from subprocess import run

ii = range(25)

for i in ii:
    run(["sbatch","-J","2wave"+str(i),"submit.cmd",str(i)])

