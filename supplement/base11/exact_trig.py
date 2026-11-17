"""Directed integer enclosures for cos(pi*r/d), used by the certificates.

Machin's identity and Taylor remainder bounds are evaluated exactly.
Every returned integer interval is in units of 2**90.
"""
from fractions import Fraction
from math import factorial, isqrt

SCALE = 1 << 90

def ceil_div(a, b):
    return -((-a)//b)

def ceil_sqrt(n):
    r = isqrt(n)
    return r + (r*r < n)

def atan_interval(q, n):
    s = sum((Fraction((-1)**j, (2*j+1)*q**(2*j+1))
             for j in range(n)), Fraction())
    t = s + Fraction((-1)**n, (2*n+1)*q**(2*n+1))
    return min(s,t), max(s,t)

def pi_interval():
    a,b = atan_interval(5,90)
    c,d = atan_interval(239,30)
    lo,hi = 16*a-4*d,16*b-4*c
    return lo.numerator*SCALE//lo.denominator, ceil_div(hi.numerator*SCALE,hi.denominator)

PI_LO, PI_HI = pi_interval()
assert 0 < PI_LO <= PI_HI and PI_HI*7 < 22*SCALE

def cos_interval(r, den):
    if not (den > 0 and abs(r) <= den):
        raise ValueError('require |r| <= den')
    r=abs(r)
    yl,yu=PI_LO*r//den,ceil_div(PI_HI*r,den)
    zl,zu=yl*yl//SCALE,ceil_div(yu*yu,SCALE)
    tl=tu=cl=cu=SCALE
    for j in range(1,21):
        div=SCALE*(2*j-1)*(2*j)
        tl,tu=tl*zl//div,ceil_div(tu*zu,div)
        if j%2:
            cl,cu=cl-tu,cu-tl
        else:
            cl,cu=cl+tl,cu+tu
    error=ceil_div(SCALE*22**42,7**42*factorial(42))
    return cl-error,cu+error
