from subprocess import run

ii = [0,2,4,5,6,7,9,10,11,12,13,14,15,16,17,19,20]

for i in ii:
    run(["sbatch","-J","2wave"+str(i),"submit.cmd",str(i)])

