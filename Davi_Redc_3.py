############# David Redction ##############

import random
import math
from sympy import O, isprime
from RNSopreration import *

p = 2 ** 250 * 3 ** 159 - 1
e = 503    
N = 22

s = 48
u = 2 * s
alpha = 0.5
l = 18
t = 10

B = [281474976710653, 281474976710651, 281474976710649, 281474976710647, 281474976710641, 281474976710639, 281474976710635,
        281474976710633, 281474976710629, 281474976710623, 281474976710621, 281474976710611, 281474976710597, 281474976710593,
        281474976710591, 281474976710581, 281474976710579, 281474976710573, 281474976710569, 281474976710567, 281474976710563, 
        281474976710561]

Q = mul(B)

def QdivBi(Q ,B):
    QdivBi = []
    for i in range(N):
        QdivBi.append(Q // B[i])
    return QdivBi

def Q1_INV(Q2, B):
    Q1_INV = []
    for i in range(N):
        Q1_INV.append(modinv(Q2[i], B[i]))
    return Q1_INV

def Q1mod_p_B(Q1, B):
    Q3 = []
    for i in range(N):
        tmp = []
        for j in range(N):
            tmp.append((Q1[i] % p) % B[j])
        Q3.append(tmp)
    return Q3

def Q1mod_p(Q1):
    Q4 = []
    for i in range(N):
        Q4.append(Q1[i] % p)
    return Q4

# applox k
def kapp(Va):
    k = int(alpha*2**l)
    for i in range(N):
        k = k + (Va[i] >> (s - l))
    k = k >> l
    return k

# r approx
def rapp(Q2, u, k):
    SS_2 = []
    sie = 0
    for i in range(N):
        sie = Q4[i]
        SS_2.append(sie>>(e - u))
    print(SS_2)
    SS_3 = multiRNS_no_red(Q2, SS_2, B, N)
    
    sum1 = 0
    for i in range(N):
        sum1=sum1+SS_3[i]
    SSq = (-((k * Q)% p)>>(e - u))
    sum1 = SSq + sum1
    quotient = (sum1<<(e-u)) //p
    return quotient

def David_reduction():
    v = []
    for i in range(N):
        v.append((ZR[i] * Q2[i]) % B[i])
    k = kapp(v)
    r = rapp(Q3, u, k)  
    vi_Qtil= []
    for i in range(N):
        vi_Qtil.append(v[i] * (Q1[i] % p))
    amodp = (sum(vi_Qtil) % p - (k * Q) % p) - r * p
    ans = []
    for i in range(N):
        ans.append(amodp % B[i])
    return ans

if __name__ == "__main__":
    X = random.randint(0, 9 * Q // 10)
    XR = INTtoRNS(X, B, N)

    Q1 = QdivBi(Q, B)
    Q2 = Q1_INV(Q1, B)
    Q3 = Q1mod_p_B(Q2, B)
    Q4 = Q1mod_p(Q1)

    c = X % p    #(CRT(z_R, B) % p)
    print('c = a * b mod p:', c)
    D_R = David_reduction()
    # print('David Reduction:', D_R)
    # print('a * b:', z_R)
    D_CRT = CRT(D_R, B)
    assert(c == D_CRT)
    print('David CRT (= c):', D_CRT)   #(CRT(z_R, B) % p)
    # print(dcrt.bit_length())
