#!/usr/bin/env python3
"""Exact certificates for the two-sided affine Hellinger contraction.

Trigonometric enclosures, square-root roundings, grid suprema and final
power comparisons use integer and Fraction arithmetic. The interpolation
constants are established in the paper's Hellinger curvature lemma.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

from exact_trig import SCALE, cos_interval, ceil_div

MESH = 200
DEN = 1 << 28
TARGETS = {
    0: Fraction(19, 20),
    1: Fraction(967, 1000),
    2: Fraction(239, 250),
    3: Fraction(481, 500),
    4: Fraction(241, 250),
}
DELTAS = {
    0: Fraction(7, 1000),
    1: Fraction(48, 10000),
    2: Fraction(6, 1000),
    3: Fraction(55, 10000),
    4: Fraction(53, 10000),
}


def ceil_sqrt(n):
    v = math.isqrt(n)
    return v if v*v == n else v+1


def probability_and_amplitude_upper(digit, mesh):
    period = 10 * mesh
    cosine_hi = []
    for t in range(period):
        r = (2*t + period) % (2*period) - period
        cosine_hi.append(cos_interval(r, period)[1])
    digits = [j for j in range(10) if j != digit]
    counts = [sum(j+h in digits for j in digits) for h in range(1, 10)]
    prob = []
    amp = []
    for t in range(period):
        numerator = 9*SCALE + sum(
            2*c*cosine_hi[h*t % period]
            for h,c in enumerate(counts, 1)
        )
        p_hi = max(0, ceil_div(numerator*DEN, 90*SCALE))
        p_hi = min(p_hi, ceil_div(9*DEN, 10))
        prob.append(p_hi)
        amp.append(ceil_sqrt(p_hi*DEN))
    return prob, amp


def profile_upper(amp, mesh):
    period = 10*mesh
    correction = ceil_div(18*DEN, 8*mesh*mesh)
    profile = [0]*period
    # rho is even. Checking one half and reflecting is exact.
    for shift in range(period//2+1):
        largest = 0
        for x in range(mesh):
            h = sum(amp[x+i*mesh] * amp[(x+i*mesh+shift) % period]
                    for i in range(10))
            largest = max(largest, h)
        value = min(DEN, ceil_div(largest, DEN)+correction)
        profile[shift] = value
        profile[-shift % period] = value
    profile[0] = DEN
    return profile


def finite_K2_max(prob, rho, mesh):
    period = 10*mesh
    rho2 = [min(DEN, ceil_div(r*r, DEN)) for r in rho]
    largest = -1
    witness = None
    # a and -a give the same supremum after kappa -> -kappa, since rho is even.
    for a in range(1, 6):
        cycle = 10//math.gcd(a,10)
        for x in range(mesh+1):
            weights = [
                sum(prob[(x+i*mesh) % period]
                    for i in range(j,10,cycle))
                for j in range(cycle)
            ]
            shifts = [(a*j*mesh) % period for j in range(cycle)]
            for kappa in range(period):
                value = sum(
                    weights[j]*rho2[(kappa+shifts[j]) % period]
                    for j in range(cycle)
                )
                if value > largest:
                    largest = value
                    witness = [a,x,kappa]
    return largest, witness


def compute(digit, mesh):
    prob, amp = probability_and_amplitude_upper(digit,mesh)
    rho = profile_upper(amp,mesh)
    num, witness = finite_K2_max(prob,rho,mesh)
    upper = Fraction(num, DEN*DEN)+Fraction(28,8*mesh*mesh)
    target = TARGETS[digit]
    delta = DELTAS[digit]
    assert upper <= target, (digit,upper,target)
    # target <= 10^(-3 delta), checked with exact integers.
    exponent = 3*delta
    n, d = exponent.numerator, exponent.denominator
    assert target.numerator**d * 10**n <= target.denominator**d
    return {
        "digit":digit, "reflected_digit":9-digit, "mesh":mesh,
        "finite_K2_numerator":num, "finite_K2_denominator":DEN*DEN,
        "grid_witness":[str(Fraction(v,mesh)) if i else v
                        for i,v in enumerate(witness)],
        "K2_upper":str(upper), "K2_target":str(target),
        "certified_delta":str(delta),
        "rho_upper_at_integers":[rho[j*mesh] for j in range(10)],
        "rho_denominator":DEN,
    }

def main():
    assert __debug__, 'Run without -O: exact assertions must remain enabled.'
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', nargs='?', type=Path,
                        default=Path(__file__).with_name('hellinger_certificates.json'))
    args = parser.parse_args()
    data = json.loads(args.path.read_text())
    assert data['status'] == 'complete'
    assert data['mesh'] == MESH and data['fixed_denominator'] == DEN
    assert [record['digit'] for record in data['digits']] == list(range(5))
    for stored in data['digits']:
        computed = compute(stored['digit'], MESH)
        assert computed == stored
        print('VERIFIED', computed['digit'], 'and', computed['reflected_digit'],
              'K2 <=', computed['K2_target'],
              'delta >=', computed['certified_delta'], flush=True)
    print('ALL FIVE HELLINGER CERTIFICATES PASSED', flush=True)

if __name__ == '__main__':
    main()
