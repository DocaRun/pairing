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

def David_reduction(z_R, Qtilinv, Qpq, Q_til_modp):
    V_R = []
    for i in range(N):
        V_R.append((z_R[i] * Qtilinv[i]) % B[i])
    k = int(alpha * l)
    r = 0
    Y = []
    for i in range (N):
        k = k + (V_R[i] >> (s - l))
        r = r + V_R[i] * (Q_til_modp[i] >> (e - u))
        Ytmp = 0
        for j in range(N):
            Ytmp = (Ytmp + V_R[i] * Qpq[i][j])
        Y.append(Ytmp % B[i])
    k = k >> l
    r = (r - (k * Q) % p) + math.floor((2 ** (e + 1)) // p) # t = 1, t is arbitrary chosen > 0
    r = r >> (1 + u)
    Z_R = []
    for i in range(N):
        rslt = (Y[i]  - (((k * Q) % p) % B[i] - (r * (p) % B[i])) % B[i])
        Z_R.append(rslt)
    return Z_R

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
    D_R = David_reduction(z_R, Qtil_inv, Qpq, Qtil_modp)
    # print('David Reduction:', D_R)
    # print('a * b:', z_R)
    dcrt =CRT(D_R, B)
    print('David CRT (= c):', dcrt)   #(CRT(z_R, B) % p)
    # print(dcrt.bit_length())
