"""Exact logarithm enclosures and affine geometry in three coordinates."""

from fractions import Fraction as Q


def log_bounds(num, den=1):
    z = Q(num-den, num+den)
    count = 100
    lo = 2*sum(z**(2*j+1)/(2*j+1) for j in range(count))
    hi = lo + 2*z**(2*count+1)/((2*count+1)*(1-z*z))
    return lo, hi


delta = Q('.096911')
kap2_range = Q('.301029'), Q('.301030')
kap5_range = 1-kap2_range[1], 1-kap2_range[0]


def add(*rows):
    return tuple(sum(entries) for entries in zip(*rows))


def scale(coef, row):
    return tuple(coef*x for x in row)


def sub(row, other):
    return add(row, scale(-1, other))


def const(x):
    return Q(x), Q(0), Q(0), Q(0)


def value(row, point):
    return sum(x*y for x, y in zip(row, (1,*point)))


def dot(x, y):
    return sum(a*b for a,b in zip(x,y))


def cross(x, y):
    a,b,c = x
    d,e,f = y
    return b*f-c*e, c*d-a*f, a*e-b*d


def intersection(rows):
    r,s,t = rows
    x,y,z = r[1:],s[1:],t[1:]
    u,v,w = cross(y,z),cross(z,x),cross(x,y)
    det = dot(x,u)
    if not det:
        return None
    return tuple(-(r[0]*a+s[0]*b+t[0]*c)/det for a,b,c in zip(u,v,w))
