import numpy as np 
import dna2mhd_utils as dn

def maxinds(fname,nmax):
    f = open(fname,'r')
    f.seek(0)
    li = []
    cs = []
    for i in range(10**5):
        a = f.readline()
        if 'Max' in a:
            b = a.split(' ')
            if 'V' in a:
                t = 'V'
            else:
                t = 'B'
            inds = [t]
            for j in b:
                c = j.split('\n')
                for k in c:
                    if k.isnumeric():
                        inds.append(int(k)-1)
            if (inds[1:] != [0,0,0]):
                li.append(inds)

    l1 = []
    cs = []
    for ik in range(len(li)):
        if li[ik] not in l1:
            l1.append(li[ik])
            cs.append(li.count(li[ik]))
        #if np.mod(ik,100) == 0: 
            # print(ik)

    lii = []
    cii = []
    for iik in range(len(l1)):
        if ((0 not in l1[iik])and(1 not in l1[iik])and(2 not in l1[iik])) and (((nmax-1) not in l1[iik]) and ((nmax-2) not in l1[iik]) and ((nmax-3) not in l1[iik])):
            lii.append(l1[iik])
            cii.append(cs[iik])
    if len(lii) > 10:
        rls = []
        a = np.argsort(cii)
        for i in a[-10:]:
            rls.append([lii[i],cii[i]])
    else:
        for i in range(len(lii)):
            rls.append([lii[i],cii[i]])

    res = []
    a = np.argsort(cs)
    for i in a[-(40-len(rls)):]:
        res.append([l1[i],cs[i]])
    res.extend(rls)
    return(res)

