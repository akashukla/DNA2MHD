import numpy as np
from numpy import format_float_positional as ff

perpscale = 0.05
parscale = 0.01

def resonancecondition(kx,ky,kz,lx,ly,lz,b1,b2,h1,h2):

    k = np.sqrt(perpscale**2.0 * (kx**2 + ky**2) + parscale**2.0 * kz**2 )
    wk = h1 * parscale * kz * (k / 2 + b1 * np.sqrt(1+(k/2)**2.0))
    l = np.sqrt(perpscale**2.0 * (lx**2 + ly**2) + parscale**2.0 * lz**2 )
    wl = h1 * parscale * lz * (l / 2 + b1 * np.sqrt(1+(l/2)**2.0))

    m = np.sqrt(perpscale**2.0 * ((kx+lx)**2 + (ky+ly)**2) + parscale**2.0 * (kz+lz)**2)
    wm = h2 * parscale * (kz+lz) * (m / 2 + b2 * np.sqrt(1+(m/2)**2.0))

    if min(np.abs(wk),np.abs(wl),np.abs(wm)) == 0:
        print(kx,ky,kz,lx,ly,lz,b1,b2,h1,h2)

    #if np.abs(wm - wk -wl) < 0.005:
    #    print(wk,wl,wm,wk+wl-wm,lx,ly,lz)
    return(np.abs(wm-wl-wk)/min(np.abs(wk),np.abs(wl),np.abs(wm)),1.0/np.abs(wm))

#kx = 0
#ky = 35
#kz = 9

kx = 0
ky = 15
kz = 9

#kx = 0
#ky = 10
#kz = 15

#kx = 0
#ky = 47
#kz = 3

#kx = 0
#ky = 40
#kz = 17

#kx = 0
#ky = 15
#kz = 19

cutoff = 63
cutoff = cutoff//3 * 2

def parallelloop(b1,b2,h1,h2):
    lx = 0
    target = 100
    best = (-99,-99,-99)
    period = 1000
    for ly in range(2,cutoff):
        for lz in range(-cutoff+1,cutoff):
            if (max(abs(kx+lx),abs(ky+ly),abs(kz+lz))< cutoff ) and (lz != 0 and lz+kz != 0) and cosineangle((kx,ky,kz),(lx,ly,lz)) < np.cos(np.pi/180):
                objective = resonancecondition(kx,ky,kz,lx,ly,lz,b1,b2,h1,h2)
                if target > objective[0]:
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
            if (max(abs(kx+lx),abs(ky+ly),abs(kz+lz))< cutoff ) and (lz != 0 and lz+kz != 0) and cosineangle((kx,ky,kz),(lx,ly,lz)) > np.cos(179*np.pi/180):
                objective = resonancecondition(kx,ky,kz,lx,ly,lz,b1,b2,h1,h2)
                if target > objective[0]:
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
                objective = resonancecondition(kx,ky,kz,lx,ly,lz,b1,b2,h1,h2)
                if target > objective[0]:
                    target = objective[0]
                    best = (lx,ly,lz)
    return(best,target,period)

def cosineangle(ks,ls):
    k = np.sqrt(perpscale**2.0 * (ks[0]**2 + ks[1]**2) + parscale**2.0 * ks[2]**2 ) 
    l = np.sqrt(perpscale**2.0 * (ls[0]**2 + ls[1]**2) + parscale**2.0 * ls[2]**2 )
    dp = perpscale**2.0 * (ks[0]*ls[0]+ks[1]*ls[1])+parscale**2.0 * ks[2]*ls[2]
    return(dp/(k*l))

def pcurleig(karr):
    pcurleig = np.cross(karr,np.array([0,0,1])).astype("complex128")
    pcurleig += 1j * np.cross(karr,pcurleig)/np.linalg.norm(karr)
    pcurleig *= 1.0/(np.sqrt(2.0) * np.linalg.norm(karr[:2]))
    return(pcurleig)

def alpha(karr):

    k = np.linalg.norm(karr)
    return(k/2 + np.sqrt(1+(k/2)**2.0))

def intkernel(karr,larr):

    # karr - interacting mode
    # larr - formed mode
    
    Ikl = np.dot(karr,pcurleig(larr-karr))*np.dot(pcurleig(karr),np.conj(pcurleig(larr)))
    # print("Int Knl ",Ikl)
    return(Ikl)

def fg(karr,larr):

    f = (1/(alpha(karr)*alpha(larr-karr)) - 1) * (np.linalg.norm(larr-karr)-np.linalg.norm(karr))/(np.linalg.norm(larr))
    g = 1/alpha(karr) - 1/alpha(larr-karr) + np.linalg.norm(larr-karr) - np.linalg.norm(karr)

    return(f,g)
    
def threetime(ks,ls):
    kx = perpscale * ks[0]
    lx = perpscale * ls[0]
    ky = perpscale * ks[1]
    ly = perpscale * ls[1]
    kz = parscale * ks[2]
    lz = parscale * ls[2]

    karr = np.array([kx,ky,kz])
    larr = np.array([lx,ly,lz])

    fh1,gh1 = fg(karr+larr,karr)
    f21,g21 = fg(larr,karr)
    fh2,gh2 = fg(karr+larr,larr)
    f12,g12 = fg(karr,larr)
    f1h,g1h = fg(karr,karr+larr)
    f2h,g2h = fg(larr,karr+larr)
    
    coeff1 = -1j/np.sqrt(1+0.25*np.linalg.norm(karr)**2.0) * ((intkernel(karr+larr,karr))*(alpha(karr)*gh1+fh1)
                                                              + (intkernel(larr,karr))*(alpha(karr)*g21+f21))/2
    
    coeff2 = -1j/np.sqrt(1+0.25*np.linalg.norm(larr)**2.0) * ((intkernel(karr+larr,larr))*(alpha(larr)*gh2+fh2)
                                                              + (intkernel(karr,larr))*(alpha(larr)*g12+f12))/2

    coeff3 = -1j/np.sqrt(1+0.25*np.linalg.norm(karr+larr)**2.0) * ((intkernel(karr,karr+larr))*(alpha(karr+larr)*g1h+f1h)
                                                                   + (intkernel(larr,karr+larr))*(alpha(karr+larr)*g2h+f2h))/2
    
    print(coeff1,coeff2,coeff3)
    threetime = 4/np.sqrt(np.abs(coeff1)*np.abs(coeff2)*np.abs(coeff3))

    return(threetime)
    
signs = [1,-1]
branches = [" Whistler "," Cyclotron "]
helicities = [" Positive "," Negative "]

for h1 in range(1):
    for h2 in range(1):
        for b1 in range(1):
            for b2 in range(1):
                best,target,period = parallelloop(signs[b1],signs[b2],signs[h1],signs[h2])
                print("Best Parallel "+helicities[h1]+branches[b1]+helicities[h2]+branches[b2],best,"\n Dot Product ",ff(cosineangle((kx,ky,kz),best),5),
                      "\n Resonance Score "+ff(target,5),"\n Period "+ff(period,5),
                      "\n Energy Time "+ff(threetime((kx,ky,kz),best),5))
                best,target,period = antiparallelloop(signs[b1],signs[b2],signs[h1],signs[h2])
                print("Best Antiparallel "+helicities[h1]+branches[b1]+helicities[h2]+branches[b2],best,"\n Dot Product ",ff(cosineangle((kx,ky,kz),best),5),
                      "\n Resonance Score "+ff(target,5),"\n Period "+ff(period,5),
                      "\n Energy Time "+ff(threetime((kx,ky,kz),best),5))
                best,target,period = perpendicularloop(signs[b1],signs[b2],signs[h1],signs[h2])
                print("Best Perpendicular "+helicities[h1]+branches[b1]+helicities[h2]+branches[b2],best,"\n Dot Product ",ff(cosineangle((kx,ky,kz),best),5),
		      "\n Resonance Score "+ff(target,5),"\n Period "+ff(period,5),
                      "\n Energy Time "+ff(threetime((kx,ky,kz),best),5))
                print("\n")
