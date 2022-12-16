# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 10:52:54 2022

@author: David1
"""

import random
import math
from RNS_operation import *
# Prime
e1 = 4
e2 = 3
e1 = 0xfa
e2 = 0x9f
p = (2**e1)*(3**e2)-1
e = 503 #size of p

# RNS Parameters and precommuted values
ma=[281474976710655,281474976710654,281474976710653,281474976710651,281474976710647,
281474976710639,281474976710633,281474976710623,281474976710617,281474976710611,
281474976710597,281474976710593,281474976710591,281474976710581,281474976710579,
281474976710573,281474976710569,281474976710567,281474976710563,281474976710561,
281474976710527,281474976710399]
Ma = ProduitRNS(ma) # compute Ma = product(ma)
n1 = len(ma)
# Ta = Ma/mai, Q1
Ta = Tai_M(Ma, ma, n1)
invTaBE = InvProdT(Ta, n1, ma) 
size = 48
u = 96

def count_bits(n):
  integer = bin(n)
  return len(integer)-2

# Fonction that computes the parameter k
# applox k
# input : Va[i] = Xa[i] * invTaBE[i] (invTaBE[i] = (Ma / ma[i]) ^ (-1))
#       : s0 = alpha = 0.1
#       : Acc = l = 18
def k_computation(Va, s0, Acc):
    k = int(s0*2**(Acc))
    for i in range(n1):
        k = k + (Va[i]>>(size-Acc))
    k = k>>Acc
    return k

# Fonction that computes Msum given V[i] = Va[i] * T^-1[i] 
# X[i] = Va[i] * (invTaBE[i] % p) % ma[j]
def RNSRed(V):
    X = []
    for j in range(n1):
        X.append(0)
        for i in range(n1):
            x = V[i] * ((Ta[i] % p) % ma[j])
            X[j] = X[j]+x
    return X

# r approx??
# input  : tab = Tai mod p (Ta = Ma / mai)
#        : e = prime bit
#        : q = u(= 96)
# output : T = tab[i] >> (s - l)
def truncTab_v2(tab, e, q):
    T = []
    n = len(tab)
    sie = 0
    for j in range(n):
        sie = tab[j]
        T.append(sie>>(e-q))  # /2**(r-q)
    return T

# r approx
# inputs : Va[i] = Xa[i] * invTaBE[i]
#        : precision = u (= 96)
#        : k
# output : q (= r) = 
def compute_quotient(Va, u, k):
    #Computes the T mod p
    Ta2 = []
    for i in range(n1):
        Ta2.append(Ta[i]%p)
    
    # Take the u most siginificant bit of T mod p
    # SS_2[i] = tab[i] >> (e - u)
    SS_2 = truncTab_v2(Ta2, e, u)

    # Does a RNS multiplication of V and SS_2, without the reduction by the ma[i]
    # SS_3[i] = Va[i] * SS_2[i]
    #         = Va[i] * tab[i] >> (e - u)
    SS_3 = multiRNS_no_red(Va, SS_2, ma, n1)
    
    # accumulation of all the term of SS_2
    # sum1[i + 1] = sum1[i] + SS_3[i]
    #         = sum1[i] + Va[i] * tab[i] >> (e - u)
    sum1 = 0
    for i in range(n1):
        sum1=sum1+SS_3[i]
    SSq = (-((k * Ma)% p)>>(e-u))
    sum1 = SSq + sum1
    # Here I use a simple // to do the division, but the "normal" way is to multiply it
    # by 8950950946455969 = 1/p * 2**(555) before a right shift.
    quotient = (sum1<<(e-u)) //p
    return quotient

# Generate a random number and convert it to RNS
# X : 元の大きな整数
# Xa: X のRNS表現
X = random.randint(0,9*Ma//10)
Xa = DecRNS(X, ma, n1)

# Computes the production V = X * T^(-1)
# inputs : Xa = X in RNS
#        : invTaBE = (Ma / mai) ^ (-1)
#        : ma = Base of RNS
#        : n1 = len(ma)
# output : Va = Xai * invTaBEi
Va = prodRNS(Xa, invTaBE, ma, n1)

# Computes Msum (here it's called SS)
# input  : Va[i] = Xa[i] * invTaBE[i]
# output : SS[i] = Va[i] * (invTaBE[i] % p) % ma[j]
SS = RNSRed(Va)

#Computes the parameter k
k = k_computation(Va, 0.1, 18)

#Computes the parameter q (= r)
q = compute_quotient(Va, u, k)

# Using k and q, we generate the 2 RNS reduction  elements noted qP and Kma
# qP[i] = (q * p) % ma[i]
qP = DecRNS(q*p, ma, n1)
# kma[i] = ((k * Ma) % p) % ma[i]
Kma =  DecRNS((k*Ma)%p, ma, n1)

# Sub Msum and Kma in RNS
# first_redutction[i] = (SS[i] - Kma[i]) % ma[i]
#                     = (Va[i] * (invTaBE[i] % p) % ma[j] - ((k * Ma) % p) % ma[i]) % ma[i]
first_reduction = sousRNS(SS,Kma,ma,n1)

# Subtract previous result to qP to get the reduced value
# out[i] = (SS[i] - Kma[i]) % ma[i] - qp[i]
out = sousRNS(first_reduction, qP, ma, n1)

#Testing part
outint = RNSDec2(out, ma, n1)
print("input :",X,"in RNS :", Xa)
print("output in integer :",outint)
print("expected  value   :",X%p)
print("correct reduction:",X % p == outint)
if (X % p == outint-p):
    print("partial reduction case")