# David Redction,a simple version

import random
import math
from sympy import O, isprime

def get_prime(x):
    while True:
        tmp = random.getrandbits(x)
        if tmp == 0: 
            continue
        if tmp == 2:
            continue
        if isprime(tmp): 
            return tmp

def xgcd(a, b):
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b != 0:
        q, a, b = a // b, b, a % b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return a, x0, y0

def modinv(a, m):
    g, x, y = xgcd(a, m)
    if g != 1:
        raise Exception('modular inverse does not exist')
    else:
        return x % m

def mul(m):
    total = 1
    for i in m:
        total *= i
    return total

def CRT(a, m):
    sum = 0
    M = mul(m)
    for i in range(len(a)):
        Mi = M // m[i]
        sum += (a[i] * modinv(Mi, m[i])) % m[i] * Mi 
    return sum % M

def multiRNS_no_red(Xa,Ya,ma,m1):
    """
    Multiplication of 2 RNS integers without the reduction by ma[i] at the end

    """
    Xr = []
    for i in range (m1):
        x = Xa[i] * Ya[i]
        Xr.append(x)
    return Xr
# Fonction that computes the parameter k
# applox k
# input : Va[i] = Xa[i] * invTaBE[i] (invTaBE[i] = (Ma / ma[i]) ^ (-1))
#       : s0 = alpha = 0.1
#       : Acc = l = 18
def k_computation(Va, s0, Acc):
    k = int(s0*2**(Acc))
    for i in range(N):
        k = k + (Va[i]>>(s-Acc))
    k = k>>Acc
    return k
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
    for i in range(N):
        Ta2.append(Qtil[i]%p)
    
    # Take the u most siginificant bit of T mod p
    # SS_2[i] = tab[i] >> (e - u)
    SS_2 = truncTab_v2(Ta2, e, u)

    # Does a RNS multiplication of V and SS_2, without the reduction by the ma[i]
    # SS_3[i] = Va[i] * SS_2[i]
    #         = Va[i] * tab[i] >> (e - u)
    SS_3 = multiRNS_no_red(Va, SS_2, B, N)
    
    # accumulation of all the term of SS_2
    # sum1[i + 1] = sum1[i] + SS_3[i]
    #         = sum1[i] + Va[i] * tab[i] >> (e - u)
    sum1 = 0
    for i in range(N):
        sum1=sum1+SS_3[i]
    SSq = (-((k * Q)% p)>>(e-u))
    sum1 = SSq + sum1
    # Here I use a simple // to do the division, but the "normal" way is to multiply it
    # by 8950950946455969 = 1/p * 2**(555) before a right shift.
    quotient = (sum1<<(e-u)) //p
    return quotient

def David_reduction(z_R, Qtil, Qtil_inv, N, B, p):
    vi = []
    for i in range(N):
        vi.append((z_R[i] * Qtil_inv[i]) % B[i])
    
    k = k_computation(vi, 0.1, 18)

    r = compute_quotient(vi, u, k)

    vi_Qtil= []
    for i in range(N):
        vi_Qtil.append(vi[i] * (Qtil[i] % p))
    amodp = (sum(vi_Qtil) % p - (k * Q) % p) - r * p

    ans = []
    for i in range(N):
        ans.append(amodp % B[i])

    return ans

if __name__ == "__main__":
    N = 22
    B = []
    s = 48
    u = 2 * s
    alpha = 0.5
    l = 18
    p2s = 2 ** s
    di1 = [3, 7, 15, 21, 27, 35, 59, 65, 77, 87, 93]
    di2 = [5, 9, 17, 23, 33, 45, 63, 75, 83, 89, 95]
    for i in range(11):
        B.append(p2s - di1[i])
        B.append(p2s - di2[i])
    Q = mul(B)
    p = 2 ** 250 * 3 ** 159 - 1
    e = 503
    a = random.getrandbits(501) %  p
    # a = 6050354745843467015840703213625286096871267259139533504340516576642729824676683714225334731411660711691098417198039283352471893134181663911193218659800
    b = random.getrandbits(501) %  p
    # b = 8810366702918147011847988616646408443159467160250681155778967775410683832587870734358528874495026302377479077378772921075066693212932226611213827176852
    a_R = []
    b_R = []
    for i in range(N):
        a_R.append(a % B[i])
        b_R.append(b % B[i])

    z_R = []
    for i in range(N):
        z_R.append((a_R[i] * b_R[i]) % B[i])

    Qtil = []
    for i in range(N):
        Qtil.append(Q // B[i])
    
    Qtil_inv = []
    for i in range(N):
        Qtil_inv.append(modinv(Qtil[i], B[i]))
    
    Qpq = []
    for i in range(N):
        Q_j = []
        for j in range(N):
            Q_j.append((Qtil[i] % p) % B[j])
        Qpq.append(Q_j)

    Qtil_modp = []
    for i in range(N):
        Qtil_modp.append(Qtil[i] % p)

    c = (a * b) % p    #(CRT(z_R, B) % p)
    print('c = a * b mod p:', c)
    D_R = David_reduction(z_R, Qtil, Qtil_inv, N, B, p)
    # print('David Reduction:', D_R)
    # print('a * b:', z_R)
    D_CRT =CRT(D_R, B)
    print('David CRT (= c):', D_CRT)   #(CRT(z_R, B) % p)
    # print(dcrt.bit_length())
    assert(c == D_CRT)
