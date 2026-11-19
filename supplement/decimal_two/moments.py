"""Shifted decimal first, fractional, and projection moment certificates."""

from fractions import Fraction as Q
from itertools import combinations

from transfer import certify, cos_interval


def check_power(base, exponent, cap):
    assert cap.numerator ** exponent.denominator < (
        cap.denominator ** exponent.denominator * base ** exponent.numerator
    )


PROFILES = {
    'moments': ((Q(1), Q('2.6789'), Q('.428')), (Q(7, 5), Q('1.8112'), Q('.258'))),
    'projection': ((Q(3, 2), Q('1.6785'), Q('.225')), (Q(4, 3), Q('1.924'), Q('.285'))),
}


def verify(*profiles):
    nodes = 6000
    grid = 10 * nodes
    cosines = [cos_interval(2 * j, grid)[1] for j in range(grid // 2 + 1)]
    representatives = [(a, b) for a, b in combinations(range(10), 2) if (a, b) <= (9 - b, 9 - a)]
    assert len(representatives) == 25
    for profile in profiles:
        for power, cap, exponent in PROFILES[profile]:
            check_power(10, exponent, cap)
            margins = [
                certify(set(range(10)) - {a, b}, cosines, power, cap)
                for a, b in representatives
            ]
            assert min(margins) > Q(1, 10**10)
            print(f'Decimal p={power}: 25 reflection representatives, cap {cap}, '
                  f'minimum rational residual {min(margins)}', flush=True)
