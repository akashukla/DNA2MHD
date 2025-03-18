import numpy as np
from numpy import format_float_positional as ff

perpscale = 0.05
parscale = 0.01

def resonancecondition(kx,ky,kz,lx,ly,lz,b1,b2,h1,h2):

    k = np.sqrt(perpscale**2.0 * (kx**2 + ky**2) + parscale**2.0 * kz**2 )
    wk = h1 * parscale * kz * (k / 2 + b1 * np.sqrt(1+k**2.0 / 4))
    l = np.sqrt(perpscale**2.0 * (lx**2 + ly**2) + parscale**2.0 * lz**2 )
    wl = h1 * parscale * lz * (l / 2 + b1 * np.sqrt(1+l**2.0 / 4))

    m = np.sqrt(perpscale**2.0 * ((kx+lx)**2 + (ky+ly)**2) + parscale**2.0 * (kz+lz)**2)
    wm = h2 * parscale * (kz+lz) * (m / 2 + b2 * np.sqrt(1+m**2.0 / 4))

    if min(np.abs(wk),np.abs(wl),np.abs(wm)) == 0:
        print(kx,ky,kz,lx,ly,lz,b1,b2,h1,h2)

    return(np.abs(wk-wl-wm)/min(np.abs(wk),np.abs(wl),np.abs(wm)),1.0/np.abs(wm))

kx = 0
ky = 35
kz = 9

#kx = 0
#ky = 15
#kz = 15

cutoff = 63
cutoff = cutoff//3 * 2

def parallelloop(b1,b2,h1,h2):
    lx = 0
    target = 100
    best = (-99,-99,-99)
    period = 1000
    for ly in range(2,cutoff):
        for lz in range(-cutoff+1,cutoff):
            if (max(abs(kx+lx),abs(ky+ly),abs(kz+lz))< cutoff ) and (lz != 0 and lz+kz != 0) and cosineangle((kx,ky,kz),(lx,ly,lz)) < 0.999:
                objective =     resonancecondition(kx,ky,kz,lx,ly,lz,b1,b2,h1,h2)
                if target >     objective[0]:
                    target = objective[0]
                    best = (lx,ly,lz)
                    period = objective[1]
    return(best,target,period)

def antiparallelloop(b1,b2,h1,h2):
    lx = 0
    target = 100
    best = (-99,-99,-99)
    for ly in range(-cutoff+1,-1):
        for lz in range(-cutoff+1,cutoff):
            if (max(abs(kx+lx),abs(ky+ly),abs(kz+lz))< cutoff ) and (lz != 0 and lz+kz != 0) and cosineangle((kx,ky,kz),(lx,ly,lz)) > -0.999:
                objective =	resonancecondition(kx,ky,kz,lx,ly,lz,b1,b2,h1,h2)
                if target >	objective[0]:
                    target = objective[0]
                    best = (lx,ly,lz)
    return(best,target,period)
    
def perpendicularloop(b1,b2,h1,h2):
    ly = 0
    target = 100
    best = (-99,-99,-99)
    for lx in range(2,cutoff-1):
        for lz in range(-cutoff+1,cutoff):
            
            if (max(abs(kx+lx),abs(ky+ly),abs(kz+lz))< cutoff ) and (lz != 0 and lz+kz != 0):
                objective =     resonancecondition(kx,ky,kz,lx,ly,lz,b1,b2,h1,h2)
                if target >     objective[0]:
                    target = objective[0]
                    best = (lx,ly,lz)
    return(best,target,period)

def cosineangle(ks,ls):
    k = np.sqrt(perpscale**2.0 * (ks[0]**2 + ks[1]**2) + parscale**2.0 * ks[2]**2 ) 
    l = np.sqrt(perpscale**2.0 * (ls[0]**2 + ls[1]**2) + parscale**2.0 * ls[2]**2 )
    dp = perpscale**2.0 * (ks[0]*ls[0]+ks[1]*ls[1])+parscale**2.0 * ks[2]*ls[2]
    return(dp/(k*l))

signs = [1,-1]
branches = [" Whistler "," Cyclotron "]
helicities = [" Positive "," Negative "]

for h1 in range(2):
    for h2 in range(2):
        for b1 in range(2):
            for b2 in range(2):
                best,target,period = parallelloop(signs[b1],signs[b2],signs[h1],signs[h2])
                print("Best Parallel "+helicities[h1]+branches[b1]+helicities[h2]+branches[b2],best,"\n Dot Product ",ff(cosineangle((kx,ky,kz),best),5),
                      "\n Resonance Score "+ff(target,5),"\n Period "+ff(period,5))
                best,target,period = antiparallelloop(signs[b1],signs[b2],signs[h1],signs[h2])
                print("Best Antiparallel "+helicities[h1]+branches[b1]+helicities[h2]+branches[b2],best,"\n Dot Product ",ff(cosineangle((kx,ky,kz),best),5),
                      "\n Resonance Score "+ff(target,5),"\n Period "+ff(period,5))
                best,target,period = perpendicularloop(signs[b1],signs[b2],signs[h1],signs[h2])
                print("Best Perpendicular "+helicities[h1]+branches[b1]+helicities[h2]+branches[b2],best,"\n Dot Product ",ff(cosineangle((kx,ky,kz),best),5),
		      "\n Resonance Score "+ff(target,5),"\n Period "+ff(period,5))
                print("\n")
