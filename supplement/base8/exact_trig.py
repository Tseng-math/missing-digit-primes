"""Directed integer enclosures for cos(pi*r/d).

Machin's identity and alternating Taylor bounds give intervals in units
of 2**90. These routines are used by the octal moment certificate.
"""
from fractions import Fraction
import math

SCALE = 1 << 90

def ceil_div(a, b):
    return -((-a) // b)

def atan_interval(q, terms):
    """Exact rational alternating-series enclosure for atan(1/q)."""
    total = sum((Fraction((-1) ** j, (2*j+1) * q ** (2*j+1))
                 for j in range(terms)), Fraction())
    next_term = Fraction((-1) ** terms, (2*terms+1) * q ** (2*terms+1))
    return min(total, total + next_term), max(total, total + next_term)

def pi_interval():
    # Machin's identity pi = 16 atan(1/5) - 4 atan(1/239).
    a, b = atan_interval(5, 90)
    c, d = atan_interval(239, 30)
    lo, hi = 16*a - 4*d, 16*b - 4*c
    return (lo.numerator * SCALE // lo.denominator,
            ceil_div(hi.numerator * SCALE, hi.denominator))

PI_LO, PI_HI = pi_interval()

def cos_interval(r, den):
    """Enclose cos(pi*r/den), |r|<=den, as fixed-point integers."""
    r = abs(r)
    yl, yu = PI_LO*r//den, ceil_div(PI_HI*r, den)
    zl, zu = yl*yl//SCALE, ceil_div(yu*yu, SCALE)
    tl = tu = SCALE
    cl = cu = SCALE
    for j in range(1, 21):
        divisor = SCALE*(2*j-1)*(2*j)
        tl, tu = tl*zl//divisor, ceil_div(tu*zu, divisor)
        if j % 2:
            cl, cu = cl-tu, cu-tl
        else:
            cl, cu = cl+tl, cu+tu
    # The Lagrange remainder after degree 40 is <= pi^42/42!.
    # pi < 22/7 is a consequence of PI_HI/SCALE < 22/7.
    assert PI_HI*7 < 22*SCALE
    remainder = ceil_div(SCALE*22**42, 7**42*math.factorial(42))
    return cl-remainder, cu+remainder
