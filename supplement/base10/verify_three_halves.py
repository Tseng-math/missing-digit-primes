#!/usr/bin/env python3
"""Exact continuous certificates for the 3/2 digital Fourier moment.

Integer trigonometric enclosures and a semiconvex interpolation bound
verify each fixed positive rational vector on the full continuous circle.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
from exact_trig import cos_interval, SCALE, ceil_div

W_DEN=1<<40
V_DEN=1<<40


def ceil_sqrt(n):
    r=math.isqrt(n)
    return r+(r*r<n)


def upper_node_weights(d,N):
    """f_d(j/(10N))^(3/2) upper bounds with denominator W_DEN."""
    M=10*N
    assert N%10==0
    # Only cos(2*pi*j/M), 0<=j<=M/2, occurs. Exact Taylor enclosures.
    cos_hi=[cos_interval(2*j,M)[1] for j in range(M//2+1)]
    digits=[j for j in range(10) if j!=d]
    counts=[sum(j+h in digits for j in digits) for h in range(1,10)]
    weights=[]
    for j in range(M):
        sq=9*SCALE
        for h,c in enumerate(counts,1):
            t=(h*j)%M
            t=min(t,M-t)
            sq+=2*c*cos_hi[t]
        sq_hi=max(0,ceil_div(sq,81))
        amp_hi=min(SCALE,ceil_sqrt(sq_hi*SCALE))
        # Round the amplitude up, then take the exact upper square root of f^3.
        aa=min(W_DEN,ceil_div(amp_hi*W_DEN,SCALE))
        ww=ceil_sqrt(ceil_div(aa**3,W_DEN))
        weights.append(min(W_DEN,ww))
    return weights




def verify_record(rec,verbose=True):
    N=rec['N'];d=rec['digit'];v=rec['vector'];g=F(rec['pressure_bound'])
    lam=F(rec['lambda'])
    assert len(v)==N and N%10==0 and min(v)>0 and max(v)<=V_DEN
    assert g>0 and lam>1
    # lambda <= 10^g, checked without evaluating logarithms.
    assert lam.numerator**g.denominator <= 10**g.numerator*lam.denominator**g.denominator
    weights=upper_node_weights(d,N)
    digest=hashlib.sha256(','.join(map(str,weights)).encode()).hexdigest()
    assert digest==rec['weight_sha256']
    maxdiff=max(abs(v[(i+1)%N]-v[i]) for i in range(N))
    # ||h||_infty<=1 and Lip(h)=N*maxdiff/V_DEN.  C=55+6 Lip(h).
    curvature_num=55*V_DEN+6*N*maxdiff
    required=lam.denominator*10*W_DEN*curvature_num
    slack_min=None
    for j in range(N):
        rem=j%10
        total=0
        for a in range(10):
            ix=a*(N//10)+j//10
            hv=(10-rem)*v[ix]+rem*v[(ix+1)%N]
            total+=weights[a*N+j]*hv
        slack=lam.numerator*v[j]*(10*W_DEN)-lam.denominator*total
        assert slack*8*N*N>=required,(d,N,j,'semiconvex slack failed')
        gap=slack*8*N*N-required
        slack_min=gap if slack_min is None else min(slack_min,gap)
    assert slack_min==int(rec['minimum_integer_margin'])
    if verbose:
        print(json.dumps({'digit':d,'N':N,'pressure_upper':str(g),
                          'lambda':str(lam),'positive_integer_margin':str(slack_min),
                          'verified':True}),flush=True)


def main():
    assert __debug__, 'Run without -O: exact assertions must remain enabled.'
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('paths', nargs='*', type=Path)
    args = parser.parse_args()
    paths = args.paths or [Path(__file__).with_name(f'moment_three_halves_d{d}.json')
                          for d in range(5)]
    for path in paths:
        verify_record(json.loads(path.read_text()))


if __name__ == '__main__':
    main()
