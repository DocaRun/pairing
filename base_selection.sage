import math
import copy

MAX = pow(2,500)

def coprime_check(T,cand):
    n = len(T)
    flag = True
    for i in range(n):
        if math.gcd(T[i],cand)!=1:
            flag = False
    return flag

def prime_factorize(n):
    if is_prime(n):
        return MAX
    fact = factor(n)
    if len(fact)>2:
        if len(fact[0])>2:
            return fact[0][0]
        else:
            return fact[1][0]
    else:
        return MAX

def eqar(A,B):
    n = 0
    if len(A)<=len(B):
        n = len(A)
    else:
        n= len(B)
    flag = True
    for i in range(n):
        if A[i]!=B[i]:
            flag = False
    return flag

def FCFSpp(w,n):
    blacklist = [0]
    previous_blacklist = [0]
    TT = []
    _round = 1
    p2w = pow(2,w)
    while True:
        _round+=1
        previous_blacklist = copy.copy(blacklist)
        i = 1
        T = []
        while len(T)<2*n:
            cand = pow(2,w)-i
            if coprime_check(T,cand) and (cand not in blacklist):
                T.append(cand)
            i+=2
        assert len(T)==2*n
        for j in range(0,2*n):
            f2 = prime_factorize(T[j])
            if f2 < (T[j]-T[2*n-1])//2:
                blacklist.append(T[j])
        TT = copy.copy(T)
        if len(previous_blacklist)==len(blacklist):
            break
            
    di = []
    for t in TT:
        di.append(p2w-t)
    di.sort()
    
    diB = []
    diC = []
    i=0
    while i<2*n:
        diB.append(di[i])
        diC.append(di[i+1])
        i+=2
    assert len(diB)==n and len(diC)==n
        
    return diB,diC

def area(w,n):
    #area_ = (25*math.ceil(w/5)*n+(w+18)*(n-1)+w*(n-1)+10*w)*n
    area_ = pow(n,2)*(180*w+126)+55*n*w
    return area_

    #return (n+10)*w+(2*w+19)*(n-1) 

if __name__ == '__main__':
    #p = 581
    #n = 10
    #w = math.ceil(p/n)
    #print('n',n)
    #print('w',w)
    #diB,diC = FCFSpp(w=w,n=n)
    #print(diB)
    #print(diC)

    p = 1006
    n = 11
    #w = math.ceil(p/n)
    w = 48
    print('n',n)
    print('w',w)
    diB,diC = FCFSpp(w,n)
    print(diB)
    print(diC)
