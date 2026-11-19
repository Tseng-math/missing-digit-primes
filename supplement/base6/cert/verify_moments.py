"""Prove uniform shifted Fourier moment bounds for all six omissions in radix 6.

Run with Python 3, with assertions enabled. All decisions use integer or
rational arithmetic. No package, input file, or saved eigenvector is needed.
"""

from fractions import Fraction as Q
from math import factorial, isqrt


NODES = 1800
RADIX = 6
CHOICES = RADIX - 1
TRIG_SCALE = 1 << 80
WEIGHT_SCALE = 1 << 32
HEIGHT = 1 << 40
PROFILES = (
    # For each moment the three entries correspond to d=0,1,2;
    # reflection covers d=5,4,3 respectively.
    (Q(1), (Q("1.98255"), Q("2.14559"), Q("2.17585")),
           (Q(".38197"), Q(".42608"), Q(".43389"))),
    (Q(4, 3), (Q("1.58803"), Q("1.66538"), Q("1.68030")),
           (Q(".25813"), Q(".28468"), Q(".28965"))),
    (Q(3, 2), (Q("1.45246"), Q("1.49930"), Q("1.50854")),
           (Q(".20833"), Q(".22604"), Q(".22947"))),
    (Q(7, 4), (Q("1.30259"), Q("1.31869"), Q("1.32199")),
           (Q(".14755"), Q(".15440"), Q(".15580"))),
)


def ceil_div(a, b):
    return -(-a // b)


def atan_interval(q, terms):
    value = sum((Q((-1)**j, (2*j+1) * q**(2*j+1))
                 for j in range(terms)), Q(0))
    adjacent = value + Q((-1)**terms, (2*terms+1) * q**(2*terms+1))
    return min(value, adjacent), max(value, adjacent)


def cosines(size=RADIX*NODES):
    # Machin's formula pi = 16 atan(1/5) - 4 atan(1/239).
    a, b = atan_interval(5, 60)
    c, d = atan_interval(239, 18)
    lo, hi = 16*a-4*d, 16*b-4*c
    scale = TRIG_SCALE
    pi_lo = (lo.numerator * scale) // lo.denominator
    pi_hi = ceil_div(hi.numerator * scale, hi.denominator)
    assert 7*pi_hi < 22*scale
    remainder = ceil_div(scale * 22**42, 7**42 * factorial(42))
    table = []
    for j in range(size//2 + 1):
        lower = 2*j*pi_lo // size
        upper = ceil_div(2*j*pi_hi, size)
        x2_lo = lower*lower // scale
        x2_hi = ceil_div(upper*upper, scale)
        term_lo = term_hi = total_hi = scale
        for k in range(1, 21):
            den = scale * (2*k-1) * (2*k)
            term_lo = term_lo * x2_lo // den
            term_hi = ceil_div(term_hi * x2_hi, den)
            total_hi += term_hi if k % 2 == 0 else -term_lo
        table.append(total_hi + remainder)
    return table


def squared_weights(digit, table, size=RADIX*NODES):
    allowed = set(range(RADIX)) - {digit}
    counts = [sum(a+j in allowed for a in allowed) for j in range(1, RADIX)]
    weights = []
    for k in range(size):
        hi = CHOICES*TRIG_SCALE
        for j, count in enumerate(counts, 1):
            r = (j*k) % size
            hi += 2*count*table[min(r, size-r)]
        w = ceil_div(hi*WEIGHT_SCALE, CHOICES**2*TRIG_SCALE)
        assert w >= 0
        weights.append(min(WEIGHT_SCALE, w))
    return weights


def upper_root(n, degree):
    """Integer ceiling root, with a dyadic shortcut."""
    if degree & (degree-1) == 0:
        v = n
        for _ in range(degree.bit_length()-1):
            v = isqrt(v)
        return v + (v**degree < n)
    lo, hi = -1, WEIGHT_SCALE
    assert hi**degree >= n
    while hi-lo > 1:
        mid = (lo+hi)//2
        if mid**degree >= n:
            hi = mid
        else:
            lo = mid
    return hi


def certify(weights, s, cap):
    a, b = s.numerator, 2*s.denominator
    scale_power = WEIGHT_SCALE**(b-a)
    powers = [upper_root(w**a * scale_power, b) for w in weights]

    def apply(v):
        # Value of the positive matrix in common units RADIX*WEIGHT_SCALE.
        out = []
        for i in range(NODES):
            total = 0
            for r in range(RADIX):
                ix = i + r*NODES
                j, rem = divmod(ix, RADIX)
                total += powers[ix] * ((RADIX-rem)*v[j] + rem*v[(j+1) % NODES])
            out.append(total)
        return out

    v = [HEIGHT]*NODES
    for _ in range(24):
        tv = apply(v)
        scale = max(tv)
        v = [HEIGHT*t//scale for t in tv]
    assert min(v) > 0 and max(v) <= HEIGHT
    lip = Q(NODES*max(abs(v[i]-v[i-1]) for i in range(NODES)), HEIGHT)
    den = RADIX*WEIGHT_SCALE
    nodal_gap = min(cap.numerator*den*y - cap.denominator*t
                    for y, t in zip(v, apply(v)))
    # Semiconvexity estimate on each node cell, proved in the accompanying proof.
    margin = Q(nodal_gap, HEIGHT*cap.denominator*den) - s*(24+4*lip)/(8*NODES**2)
    assert margin > 0, margin
    return margin


def check_exponent(cap, exponent):
    a, b = exponent.numerator, exponent.denominator
    assert cap.numerator**b < RADIX**a * cap.denominator**b


def main():
    if not __debug__:
        raise RuntimeError("Run this verifier with assertions enabled.")
    assert NODES % (2*RADIX) == 0
    # Outward decimal enclosures used for delta and the costs of the
    # two pure-prime exponents in the later comparisons.
    denominator = 10**6
    assert 6**898244 < 5**denominator < 6**898245
    assert 6**386852 < 2**denominator < 6**386853
    table = cosines()
    for d in range(RADIX//2):
        weights = squared_weights(d, table)
        for s, caps, exponents in PROFILES:
            cap, exponent = caps[d], exponents[d]
            margin = certify(weights, s, cap)
            check_exponent(cap, exponent)
            assert margin > Q(1, 100000)
            print(f"d={d} or {5-d}, p={s}, lambda={cap}, exponent={exponent}, "
                  f"minimum residual={margin}")


if __name__ == "__main__":
    main()
