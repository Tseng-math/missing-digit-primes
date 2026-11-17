#!/usr/bin/env python3
"""All six base-11 Hellinger profiles, with directed integer arithmetic.

No stored profile, floating point, numerical linear algebra, or external
package is needed. The interpolation corrections 18/(8M^2) and 20/(8M^2)
are proved in the paper. Each calculation covers all real states.
"""
from fractions import Fraction as F
from exact_trig import SCALE, ceil_div, ceil_sqrt, cos_interval

BASE=11
MESH=80
DEN=1<<28
TARGET=F(24,25)
DELTA=F(1,250)

def probability_upper(digit):
    period=BASE*MESH
    cs=[]
    for t in range(period):
        r=(2*t+period)%(2*period)-period
        cs.append(cos_interval(r,period)[1])
    ds=[a for a in range(BASE) if a != digit]
    counts=[sum(a+h in ds for a in ds) for h in range(1,BASE)]
    prob=[];amp=[]
    for t in range(period):
        num=10*SCALE+sum(2*c*cs[h*t%period] for h,c in enumerate(counts,1))
        v=max(0,min(ceil_div(10*DEN,11),ceil_div(num*DEN,110*SCALE)))
        prob.append(v);amp.append(ceil_sqrt(v*DEN))
    return prob,amp

def profile_upper(amp):
    period=BASE*MESH
    correction=ceil_div(18*DEN,8*MESH*MESH)
    profile=[0]*period
    for shift in range(period//2+1):
        largest=0
        for x in range(MESH):
            h=sum(amp[x+i*MESH]*amp[(x+i*MESH+shift)%period] for i in range(BASE))
            largest=max(largest,h)
        value=min(DEN,ceil_div(largest,DEN)+correction)
        profile[shift]=profile[-shift%period]=value
    profile[0]=DEN
    return profile

def check(digit):
    prob,amp=probability_upper(digit)
    rho=profile_upper(amp)
    rho2=[min(DEN,ceil_div(r*r,DEN)) for r in rho]
    period=BASE*MESH
    largest=-1;witness=None
    for a in range(1,BASE//2+1):
        for x in range(MESH+1):
            ps=[prob[(x+i*MESH)%period] for i in range(BASE)]
            shifts=[a*i*MESH%period for i in range(BASE)]
            for offset in range(period):
                value=sum(ps[i]*rho2[(offset+shifts[i])%period] for i in range(BASE))
                if value>largest:
                    largest=value;witness=(a,x,offset)
    upper=F(largest,DEN*DEN)+F(20,8*MESH*MESH)
    assert upper<=TARGET,(digit,upper,TARGET)
    power=3*DELTA
    assert TARGET.numerator**power.denominator*BASE**power.numerator<=TARGET.denominator**power.denominator
    print('VERIFIED base 11 digits',digit,10-digit,'K2 upper',upper,'<',TARGET,
          'delta',DELTA,'witness',witness,flush=True)
    return upper

def main():
    assert __debug__, 'Run without -O: assertions are proof checks.'
    # Parseval derivative estimates for A_i(x)=110^(-1/2) P_d(e((i+x)/11))
    # after multiplication by the centering phase e(-5(i+x)/11).
    assert (2*F(22,7)/11)**2*110/10 < 4
    assert (2*F(22,7)/11)**4*1958/10 < 25
    for digit in range(6):
        check(digit)
    print('ALL SIX BASE-11 HELLINGER CERTIFICATES PASSED',flush=True)

if __name__=='__main__':
    main()
