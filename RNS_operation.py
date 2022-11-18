# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 10:02:34 2021

@author: David1
"""
import math

def ProduitRNS (ma):
    """
    Computes Ma = products of all ma(i)

    Parameters
    ----------
    ma : RNS base

    Returns
    -------
    Ma : Integer
    Ma = products(ma(i))

    """
    m1 = len(ma)
    Ma = 1
    for i in range (m1):
        Ma = ma[i] * Ma
    return Ma

def Tai_M (Ma,ma,m1):
    """
    Computes (Ta) = (Ma/ma(i))

    Parameters
    ----------
    Ma : Integer
    Product integer Ma
    ma : TYPE
    RNS base
    m1 : integer
    Size of RNS base ma

    Returns
    -------
    Ta : List

    """
    Ta = []
    for i in range (m1):
        Ta.append(Ma//ma[i])
    return Ta

def egcd(a, b):
    """
    Extended Euclidean Algorithm

    Parameters
    ----------
    a : integer
    
    b : integer

    """
    if a == 0:
        return (b, 0, 1)
    else :
        gcd, x, y = egcd(b % a, a)
        return (gcd, y - (b//a) * x, x)

def DecRNS (X,ma,m1):
    """
    Converter from Decimal to RNS representation

    Parameters
    ----------
    X : integer
    
    ma : RNS base
    
    m1 : integer


    Returns
    -------
    Xa : integer
    X in the RNS representation using the base ma

    """
    Xa = []
    for i in range (m1):
        Xa.append((X%ma[i]))
    return Xa

def RNSDec (Xd,md,m4,Td,Md):
    """
    Converter from RNS representation to Decimale

    Parameters
    ----------
    Xd : integer
    Integer to convert
    md : RNS Base
    
    m4 : integer
    len(md)
    Td : list
    Td = (Md/md(i))
    Md : integer
    product of all elements in the RNS base

    Returns
    -------
    X : integer
    Integer xd in decimale.

    """
    X = 0
    for i in range (m4):
        pgcd, n, inv = egcd (Td[i],md[i])
        x = (Xd[i] * n) % md[i]
        X = X + (x * Td[i]) 
    X = X % Md
    return X

def RNSDec2 (Xd,md,m4):
    """
    Converter from RNS representation to Decimale with only 3 parameters

    """
    X = 0
    Md = ProduitRNS(md)
    Td = Tai_M(Md, md, m4)
    for i in range (m4):
        pgcd, n, inv = egcd (Td[i],md[i])
        x = (Xd[i] * n) % md[i]
        X = X + (x * Td[i]) 
    X = X % Md
    return X

def InvProdT(Ta,m1,ma):
    """
    Compute the RNS of T^-1    
    
    """
    invT = []
    for i in range (m1):
        pgcd, inv, n = egcd (Ta[i],ma[i])
        invT.append(inv % ma[i])
    return invT


#Arithmetic operants functions

def addRNS(Xa,Ya,ma,m1):
    """
    Addition of 2 RNS integers

    """
    Xr = []
    #print(m1)
    for i in range (m1):
        x = Xa[i] + Ya[i]
        Xr.append(x% ma[i])
    return Xr

def sousRNS(Xa,Ya,ma,m1):
    """
    Subtraction of 2 RNS integers

    """
    Xr = []
    #print(Xa,Ya,ma,m1)
    for i in range (m1):
        # if Xa[i]-Ya[i]>0:
        #     x = Xa[i] - Ya[i]
        # else:
        #     x = Xa[i] - Ya[i] + ma[i]
        x = Xa[i] - Ya[i]
        Xr.append(x% ma[i])
    #print(Xr)
    return Xr

def prodRNS(Xa,Ya,ma,m1):
    """
    Multiplication of 2 RNS integers

    """
    Xr = []
    for i in range (m1):
        x = Xa[i] * Ya[i]
        Xr.append(x % ma[i])
    return Xr
def multiRNS_no_red(Xa,Ya,ma,m1):
    """
    Multiplication of 2 RNS integers without the reduction by ma[i] at the end

    """
    Xr = []
    for i in range (m1):
        x = Xa[i] * Ya[i]
        Xr.append(x)
    return Xr
