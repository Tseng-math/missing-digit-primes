"""Verify the continuous octal first-moment supersolutions."""

from fractions import Fraction as F
from math import factorial, isqrt
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / 'first_moment_nodes.json'
S = 10**24
LAM = F(2769, 10000)

def ceildiv(a,b):
    return -((-a)//b)

def rational(x):
    return f'{x.numerator}/{x.denominator}'

def atan_interval(q,terms):
    v=sum((F((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(terms)),F(0))
    w=F((-1)**terms,(2*terms+1)*q**(2*terms+1))
    return min(v,v+w),max(v,v+w)

def pi_interval():
    a,b=atan_interval(5,50)
    c,d=atan_interval(239,15)
    lo,hi=16*a-4*d,16*b-4*c
    assert F(3)<lo<hi<F(22,7)
    ilo=lo.numerator*S//lo.denominator
    ihi=ceildiv(hi.numerator*S,hi.denominator)
    assert ihi-ilo<=1
    return ilo,ihi

def cosine_upper_table(den):
    p0,p1=pi_interval()
    terms=24
    rem=ceildiv(S*4**(2*terms+2),factorial(2*terms+2))
    table=[]
    for k in range(den//2+1):
        x0=p0*(2*k)//den
        x1=ceildiv(p1*(2*k),den)
        y0=x0*x0//S
        y1=ceildiv(x1*x1,S)
        a0=a1=c0=c1=S
        for j in range(1,terms+1):
            div=S*(2*j-1)*(2*j)
            a0=a0*y0//div
            a1=ceildiv(a1*y1,div)
            if j%2:
                c0-=a1;c1-=a0
            else:
                c0+=a0;c1+=a1
        assert c0-rem<=c1+rem
        table.append(min(S,c1+rem))
    return table

def verify_one(row,coshi):
    assert row['base']==8
    d=row['forbidden_digit']
    assert d in range(4)
    n0,scale0=row['coarse_grid_size'],row['denominator']
    ref,coarse=row['refinement'],row['numerators']
    assert n0>=2 and scale0>0 and ref>=1 and len(coarse)==n0
    assert all(isinstance(v,int) and v>0 for v in coarse)
    g=[(ref-r)*coarse[i]+r*coarse[(i+1)%n0] for i in range(n0) for r in range(ref)]
    n,scale=n0*ref,scale0*ref
    den=8*n
    assert len(coshi)==den//2+1
    counts=[8-r-int(d+r<=7)-int(d-r>=0) for r in range(1,8)]
    assert sum(counts)==21 and all(c>=0 for c in counts)
    worst=None;wi=None
    residual_den=64*S*scale*LAM.denominator
    for i in range(n):
        total=0
        for j in range(8):
            ell=i+j*n
            square=7*S
            for r,count in enumerate(counts,1):
                k=r*ell%den
                k=min(k,den-k)
                square+=2*count*coshi[k]
            square_hi=min(S,max(0,ceildiv(square,49)))
            target=square_hi*S
            w=isqrt(target)
            if w*w<target:w+=1
            assert w*w>=target
            ix,rr=divmod(ell,8)
            interpolated=(8-rr)*g[ix]+rr*g[(ix+1)%n]
            total+=w*interpolated
        rn=total*LAM.denominator-64*S*LAM.numerator*g[i]
        if worst is None or rn>worst:worst,wi=rn,i
    residual=F(worst,residual_den)
    minimum,maximum=F(min(g),scale),F(max(g),scale)
    lg=F(n*max(abs(g[(i+1)%n]-g[i]) for i in range(n)),scale)
    # Removing the common phase centred at 7/2 gives sum |a-7/2|<=16.
    # Hence Lip(f_d)<=32*pi/7<704/49, uniformly in d.
    lf=F(704,49)
    lr=(lf*maximum+lg)/8+LAM*lg
    guard=lr/(2*n)
    margin=-residual-guard
    assert margin>0,(d,rational(margin))
    return dict(digit=d,reflected_digit=7-d,grid=n,
        min_g=rational(minimum),max_g=rational(maximum),lip_g=rational(lg),
        maximum_node_residual=rational(residual),maximum_node=wi,
        continuous_guard=rational(guard),strict_margin=rational(margin))


def main():
    rows = json.loads(INPUT.read_text())['rows']
    assert len(rows) == 4 and {row['forbidden_digit'] for row in rows} == set(range(4))
    assert all(row['coarse_grid_size'] == 4096 and row['refinement'] == 4 for row in rows)
    assert 2769**400 * 8**247 < 10000**400
    cosine_bounds = cosine_upper_table(8 * 16384)
    for row in rows:
        result = verify_one(row, cosine_bounds)
        print(f"Digit {result['digit']}: residual margin {result['strict_margin']}", flush=True)
    print('Octal first-moment certificates verified for all eight digits.')


if __name__ == '__main__':
    main()
