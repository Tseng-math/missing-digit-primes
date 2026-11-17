#!/usr/bin/env python3
"""Exact continuous-state base-eleven Fourier moment certificates.

Each stored rational supersolution is verified on the full circle using
directed trigonometric enclosures and the semiconvex interpolation lemma.
"""
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path
from exact_trig import SCALE, ceil_div, ceil_sqrt, cos_interval

BASE=11
DEN=1 << 40
NODES=1100
SPECS={
    'one':(F(1),F(27,77),F(23,10),41,4),
    'fractional':(F(235,154),F(59,433),F(137,100),62,6),
    'three_halves':(F(3,2),F(287,2000),F(7,5),62,6),
}
HERE=Path(__file__).resolve().parent

@lru_cache(None)
def cosine_table(N):
    M=BASE*N
    return tuple(cos_interval(2*j,M)[1] for j in range(M//2+1))

@lru_cache(None)
def amplitudes(digit,N):
    M=BASE*N
    cs=cosine_table(N)
    ds=[a for a in range(BASE) if a != digit]
    counts=[sum(a+h in ds for a in ds) for h in range(1,BASE)]
    answer=[]
    for j in range(M//2+1):
        sq=10*SCALE
        for h,c in enumerate(counts,1):
            t=h*j%M
            sq += 2*c*cs[min(t,M-t)]
        sq=max(0,ceil_div(sq,100))
        a=min(SCALE,ceil_sqrt(sq*SCALE))
        answer.append(min(DEN,ceil_div(a*DEN,SCALE)))
    return tuple(answer)

def power_ok(a,w,p):
    return (0 <= a <= DEN and 0 <= w <= DEN
            and w**p.denominator*DEN**(p.numerator-p.denominator)>=a**p.numerator)

def full_weights(half,N):
    M=BASE*N
    assert len(half)==M//2+1
    return [half[min(j,M-j)] for j in range(M)]

def margin(weights,v,N,lam,c0,c1):
    md=max(abs(v[(j+1)%N]-v[j]) for j in range(N))
    curvature=c0*DEN+c1*N*md
    required=lam.denominator*BASE*DEN*curvature
    smallest=None
    for j in range(N):
        rem=j%BASE
        total=0
        for a in range(BASE):
            ix=a*(N//BASE)+j//BASE
            hv=(BASE-rem)*v[ix]+rem*v[(ix+1)%N]
            total += weights[a*N+j]*hv
        slack=lam.numerator*v[j]*BASE*DEN-lam.denominator*total
        gap=slack*8*N*N-required
        smallest=gap if smallest is None else min(smallest,gap)
    return smallest

def verify_record(r):
    kind=r['kind']; d=r['digit']; N=r['N']
    assert kind in SPECS and type(d) is int and 0 <= d <= 5
    assert N==NODES and N%BASE==0
    p,g,lam,c0,c1=SPECS[kind]
    assert r['p']==str(p) and r['g']==str(g) and r['lambda']==str(lam)
    assert r['base']==BASE and r['denominator']==DEN
    assert r['curvature']==[c0,c1]
    # Centering the exponents gives the two rational derivative bounds below.
    # The slightly larger common bounds are shared with the decimal proof.
    assert F(2)*F(22,7)*30/(10*11) <= F(110,63)
    assert (2*F(22,7))**2*110/(10*11**2) <= F(2662,735)
    assert c0 >= BASE*p*F(2662,735) and c1 >= 2*p*F(110,63)
    assert lam.numerator**g.denominator <= BASE**g.numerator*lam.denominator**g.denominator
    v=r['vector']; half=r['half_weights']
    assert len(v)==N and all(type(x) is int and 0 < x <= DEN for x in v)
    amps=amplitudes(d,N)
    assert len(half)==len(amps)
    assert all(type(w) is int and power_ok(a,w,p) for a,w in zip(amps,half))
    gap=margin(full_weights(half,N),v,N,lam,c0,c1)
    assert gap>0 and str(gap)==r['minimum_integer_margin']
    print('VERIFIED base 11',kind,'digits',d,BASE-1-d,'margin',gap,flush=True)
    return gap

def main():
    assert __debug__, 'Run without -O: exact assertions must remain enabled.'
    records = json.loads((HERE/'moments.json').read_text())
    assert [(r['kind'],r['digit']) for r in records] == [
        (kind,d) for kind in SPECS for d in range(6)]
    for record in records:
        verify_record(record)
    print('ALL 18 BASE-11 CONTINUOUS MOMENT CERTIFICATES PASSED', flush=True)

if __name__ == '__main__':
    main()
