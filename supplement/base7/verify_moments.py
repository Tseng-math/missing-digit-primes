"""Prove uniform shifted Fourier moment bounds for all seven omissions.

Run with Python 3, with assertions enabled. All decisions use integer or
rational arithmetic. No package, input file, or saved eigenvector is needed.
"""

from fractions import Fraction as Q
from math import factorial, isqrt


NODES = 2100
TRIG_SCALE = 1 << 80
WEIGHT_SCALE = 1 << 32
HEIGHT = 1 << 40
PROFILES = (
    # (moment, upper eigenvalue, exponent on the full shifted grid)
    (Q(1), Q("2.17457"), Q(".39922")),
    (Q(4, 3), Q("1.64337"), Q(".2553")),
    (Q(7, 5), Q("1.56675"), Q(".23075")),
    (Q(3, 2), Q("1.46669"), Q(".19684")),
    # Last exponent is interpreted relative to the exact Parseval exponent.
    (Q(127, 64), Q("1.17229"), Q(".159") / 64),
)


def ceil_div(a, b):
    return -(-a // b)


def atan_interval(q, terms):
    value = sum((Q((-1)**j, (2*j+1) * q**(2*j+1))
                 for j in range(terms)), Q(0))
    adjacent = value + Q((-1)**terms, (2*terms+1) * q**(2*terms+1))
    return min(value, adjacent), max(value, adjacent)


def cosines():
    # Machin's formula pi = 16 atan(1/5) - 4 atan(1/239).
    a, b = atan_interval(5, 60)
    c, d = atan_interval(239, 18)
    lo, hi = 16*a-4*d, 16*b-4*c
    scale = TRIG_SCALE
    pi_lo = (lo.numerator * scale) // lo.denominator
    pi_hi = ceil_div(hi.numerator * scale, hi.denominator)
    assert 7*pi_hi < 22*scale
    size = 7*NODES
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


def squared_weights(digit, table):
    size = 7*NODES
    allowed = set(range(7)) - {digit}
    counts = [sum(a+j in allowed for a in allowed) for j in range(1, 7)]
    weights = []
    for k in range(size):
        hi = 6*TRIG_SCALE
        for j, count in enumerate(counts, 1):
            r = (j*k) % size
            hi += 2*count*table[min(r, size-r)]
        w = ceil_div(hi*WEIGHT_SCALE, 36*TRIG_SCALE)
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
        # Value of the positive matrix in common units 7*WEIGHT_SCALE.
        out = []
        for i in range(NODES):
            total = 0
            for r in range(7):
                ix = i + r*NODES
                j, rem = divmod(ix, 7)
                total += powers[ix] * ((7-rem)*v[j] + rem*v[(j+1) % NODES])
            out.append(total)
        return out

    v = [HEIGHT]*NODES
    for _ in range(24):
        tv = apply(v)
        scale = max(tv)
        v = [HEIGHT*t//scale for t in tv]
    assert min(v) > 0 and max(v) <= HEIGHT
    lip = Q(NODES*max(abs(v[i]-v[i-1]) for i in range(NODES)), HEIGHT)
    den = 7*WEIGHT_SCALE
    nodal_gap = min(cap.numerator*den*y - cap.denominator*t
                    for y, t in zip(v, apply(v)))
    # Semiconvexity estimate on each node cell, proved in the integrated section.
    margin = Q(nodal_gap, HEIGHT*cap.denominator*den) - s*(28+4*lip)/(8*NODES**2)
    assert margin > 0, margin
    return margin


def check_exponent(cap, exponent):
    a, b = exponent.numerator, exponent.denominator
    assert cap.numerator**b < 7**a * cap.denominator**b


def certify_third(weights):
    # Here h=1 suffices. If w/WEIGHT_SCALE is an upper square amplitude,
    # y^2*WEIGHT_SCALE >= w^3 gives an upper cubic amplitude y/WEIGHT_SCALE.
    powers = [upper_root(ceil_div(w**3, WEIGHT_SCALE), 2) for w in weights]
    maximum = max(sum(powers[i + r*NODES] for r in range(7))
                  for i in range(NODES))
    cap = Q(103, 100)
    margin = cap-Q(maximum, WEIGHT_SCALE)-Q(3*28, 8*NODES**2)
    assert margin > 0
    check_exponent(cap, Q(".016"))
    return margin


def main():
    if not __debug__:
        raise RuntimeError("Run this verifier with assertions enabled.")
    assert NODES % 14 == 0
    table = cosines()
    scores = [[] for _ in PROFILES]
    third = []
    for d in range(4):
        weights = squared_weights(d, table)
        third.append(certify_third(weights))
        for i, (s, cap, exponent) in enumerate(PROFILES):
            scores[i].append(certify(weights, s, cap))
    for (s, cap, exponent), margins in zip(PROFILES, scores):
        # 7^delta = 7/6 in the final row.
        check_exponent(cap * (Q(6, 7) if s == Q(127, 64) else 1), exponent)
        assert min(margins) > Q(26, 10000000)
        if s in (Q(1), Q(4, 3), Q(3, 2)):
            assert min(margins) > Q(1, 100000)
        print(f"p={s}, lambda={cap}, exponent allowance={exponent}, "
              f"minimum residual={min(margins)}")
    print(f"p=3, lambda=103/100, exponent allowance=2/125, "
          f"minimum residual={min(third)}")


if __name__ == "__main__":
    main()
