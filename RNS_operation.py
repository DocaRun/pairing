# -*- coding: utf-8 -*-
"""
create : 2022/12/1
@author: DocaRun
"""
import math

def RNS (X, B):
    Xa = []
    for i in B:
        Xa.append(X % i)
    return Xa

def egcd(a, b):
    """
    Extended Euclidean Algorithm

    """
    if a == 0:
        return (b, 0, 1)
    else :
        gcd, x, y = egcd(b % a, a)
        return (gcd, y - (b//a) * x, x)
