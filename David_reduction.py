"""
create : 2022/12/1
@author: DocaRun
"""

import random
import math

# Prime
e1 = 0xfa
e2 = 0x9f
p = (2 ** e1) * (3 ** e2) - 1
e = 503

# RNS Parameters and precommuted values
B = [281474976710655,281474976710654,281474976710653,281474976710651,281474976710647,
     281474976710639,281474976710633,281474976710623,281474976710617,281474976710611,
     281474976710597,281474976710593,281474976710591,281474976710581,281474976710579,
     281474976710573,281474976710569,281474976710567,281474976710563,281474976710561,
     281474976710527,281474976710399]
Q = math.prod(B)
N = len(B)

alpha = 0.1
s = 48
l = 18
u = 96

def RNS (X, B):
    Xa = []
    for i in B:
        Xa.append(X % i)
    return Xa

def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else :
        gcd, x, y = egcd(b % a, a)
        return (gcd, y - (b//a) * x, x)

def kapp(Va):
    k = int(alpha * (2 ** l))
    for i in range(N):
        k = k + (Va[i] >> (s - l))
    k = k >> l
    return k

def rapp(Va, u, k, Q1):
    r = 0
    for i in range(N):
        r += Va[i] * ((Q1[i] % p) >> (e - u))
    r -= ((k * Q) % p) >> (e - u)
    # Here I use a simple // to do the division, but the "normal" way is to multiply it
    # by 8950950946455969 = 1/p * 2**(555) before a right shift. 
    # @David
    r = (r << (e - u)) // p
    return r

def mainReduction():
    X = random.randint(0, 9 * Q // 10)
    XR = RNS(X, B)
    Q1 = []
    for i in B:
        Q1.append(Q // i)
    Q2 = []
    for i in range (N):
        pgcd, inv, n = egcd (Q1[i], B[i])
        Q2.append(inv % B[i])
    Va = []
    for i in range (N):
        Va.append((XR[i] * Q2[i]) % B[i])
    Z = []
    for j in range(N):
        Z.append(0)
        for i in range(N):
            Z[j] = Z[j] + Va[i] * ((Q1[i] % p) % B[j])
    # k,r approximation
    k = kapp(Va)
    r = rapp(Va, u, k, Q1)
    kQR = RNS((k * Q) % p, B)
    rpR = RNS(r * p, B)

    out = []
    for i in range (N):
        tmp = Z[i] - kQR[i] - rpR[i]
        out.append(tmp % B[i])
    outint = 0
    for i in range (N):
        pgcd, n, inv = egcd (Q1[i],B[i])
        x = (out[i] * n) % B[i]
        outint = outint + (x * Q1[i]) 
    outint = outint % Q
    return X, XR, outint


if __name__ == "__main__":
    X, XR, outint = mainReduction()
    print("input :", X, "in RNS :", XR)
    print("output in integer :", outint)
    print("expected  value   :", X % p)
    print("correct reduction :", X % p == outint)
    if (X % p == outint-p):
        print("partial reduction case")
    
    ###### test ######
    succ = 0
    pati = 0
    for i in range(10000):
        X, XR, outint = mainReduction()
        # print('X',X % p)
        # print('XR', XR)
        # print('outint',outint)
        # print('i', i)
        if (X % p == outint):
            succ += 1
        elif(X % p == outint - p):
            pati += 1
        else:
            exit()
    print('The rate of correct reduction case: {} %'.format(succ / 10000 * 100))
    print('The rate of partial reduction case: {} %'.format(pati / 10000 * 100))
    print('The rate of uncorrect reduction case: {} %'.format((10000 - succ - pati) / 10000 * 100))