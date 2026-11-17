#!/usr/bin/env python3
"""Exact continuous certificates for the classical digital moment bounds.

The moment orders are 1 and 235/154. The stored positive rational vectors
are verified using integer trigonometric enclosures and the semiconvex
interpolation bound in the paper.
"""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import argparse
import hashlib
import json
import math

from exact_trig import SCALE, cos_interval, ceil_div

DEN = 1 << 40
NODES = 1000
SPECS = {
    "one": (F(1), F(27, 77), F(1121, 500), 37, 4),
    "classical": (F(235, 154), F(59, 433), F(1368543, 1000000), 56, 6),
}


def ceil_sqrt(n):
    s = math.isqrt(n)
    return s + (s*s < n)


@lru_cache(maxsize=None)
def cosine_upper_table(N):
    M = 10*N
    return tuple(cos_interval(2*j, M)[1] for j in range(M//2+1))


def amplitude_upper_nodes(digit, N):
    M = 10*N
    cosines = cosine_upper_table(N)
    digits = [j for j in range(10) if j != digit]
    counts = [sum(j+h in digits for j in digits) for h in range(1, 10)]
    amplitudes = []
    for j in range(M//2+1):
        upper_squared = 9*SCALE
        for h,c in enumerate(counts, 1):
            t = h*j % M
            upper_squared += 2*c*cosines[min(t,M-t)]
        upper_squared = max(0, ceil_div(upper_squared,81))
        upper_amplitude = min(SCALE,ceil_sqrt(upper_squared*SCALE))
        amplitudes.append(min(DEN,ceil_div(upper_amplitude*DEN,SCALE)))
    return amplitudes


def verify_weight_upper(amplitude, weight, p):
    assert 0 <= amplitude <= DEN and 0 <= weight <= DEN
    # weight / DEN >= (amplitude / DEN)^p.
    assert p >= 1
    return weight**p.denominator * DEN**(p.numerator-p.denominator) >= amplitude**p.numerator


def full_weights(half, N):
    M = 10*N
    assert len(half) == M//2+1
    return [half[min(j,M-j)] for j in range(M)]


def minimum_margin(weights, vector, N, lam, c0, c1):
    maxdiff = max(abs(vector[(j+1)%N]-vector[j]) for j in range(N))
    curvature_num = c0*DEN + c1*N*maxdiff
    required = lam.denominator*10*DEN*curvature_num
    margin = None
    for j in range(N):
        rem = j%10
        total = 0
        for a in range(10):
            ix = a*(N//10)+j//10
            interpolated = (10-rem)*vector[ix]+rem*vector[(ix+1)%N]
            total += weights[a*N+j]*interpolated
        slack = lam.numerator*vector[j]*10*DEN-lam.denominator*total
        gap = slack*8*N*N-required
        margin = gap if margin is None else min(margin,gap)
    return margin


def verify_record(record, verbose=True):
    kind = record['kind']
    assert kind in SPECS
    p,g,lam,c0,c1 = SPECS[kind]
    d = record['digit']
    N = record['N']
    assert type(d) is int and 0 <= d <= 4
    assert N == NODES and N%10 == 0
    assert record['p'] == str(p)
    assert record['pressure_bound'] == str(g)
    assert record['lambda'] == str(lam)
    assert record['denominator'] == DEN
    assert record['curvature_constants'] == [c0,c1]
    assert c0 >= 10*p*F(2662,735)
    assert c1 >= 2*p*F(110,63)
    assert lam.numerator**g.denominator <= 10**g.numerator*lam.denominator**g.denominator
    v = record['vector']
    half = record['half_weights']
    assert len(v) == N and all(type(x) is int and 0 < x <= DEN for x in v)
    amp = amplitude_upper_nodes(d,N)
    assert len(half) == len(amp)
    assert all(type(w) is int and verify_weight_upper(a,w,p) for a,w in zip(amp,half))
    weights = full_weights(half,N)
    margin = minimum_margin(weights,v,N,lam,c0,c1)
    assert margin > 0,(kind,d,'nonpositive continuous-state margin',margin)
    assert str(margin) == record['minimum_integer_margin']
    assert record['status'] == 'exact_upper_verified'
    if verbose:
        print('VERIFIED',kind,'digits',d,'and',9-d,'p =',p,'g =',g,
              'positive_integer_margin =',margin,flush=True)
    return {'kind':kind,'digit':d,'reflected_digit':9-d,
            'p':str(p),'g':str(g),'lambda':str(lam),
            'positive_integer_margin':str(margin),'verified':True}




def verify_all(directory=None):
    assert __debug__, 'Do not disable exact checks with -O.'
    directory = Path(directory) if directory else Path(__file__).resolve().parent
    results = []
    for kind in SPECS:
        for d in range(5):
            path = directory/f'moment_{kind}_d{d}.json'
            result = verify_record(json.loads(path.read_text()))
            result['certificate_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
            results.append(result)
    print('ALL 10 CLASSICAL MOMENT CERTIFICATES PASSED',flush=True)
    return results

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    verify_all(args.directory)

if __name__ == '__main__':
    main()
