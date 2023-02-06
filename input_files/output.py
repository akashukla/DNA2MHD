fname = 'DNAHD.out10884473'

import numpy as np 

def maxinds(fname):

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
    print(len(li)/8)
    l1 = []
    cs = []
    for ik in range(len(li)):
        if li[ik] not in l1:
            l1.append(li[ik])
            cs.append(li.count(li[ik]))
        if np.mod(ik,100) == 0: 
           print(ik)
    res = l1
    i = 1
    while (len(res) > 40):
        print(i)
        res = []
        for k in range(len(l1)):
            if (cs[len(l1)-1-k] >= i) and ((l1[len(l1)-1-k],cs[len(l1)-1-k]) not in res):
                res.append((l1[len(l1)-1-k],cs[len(l1)-1-k]))
        i += 1
    print(len(res))
    for ind in res:
        print(ind)
    return(res)

dum = maxinds(fname)
