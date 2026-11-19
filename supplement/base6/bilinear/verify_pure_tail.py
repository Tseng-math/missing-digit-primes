"""Directed bounds for the six-position joint coset mean of the digit transform."""
from runpy import run_path
from pathlib import Path
from fractions import Fraction as Q
from math import isqrt

math_tools = run_path(str(Path(__file__).resolve().parents[1] /
                          "cert/verify_moments.py"))
SCALE = math_tools["TRIG_SCALE"]
ceil_div = math_tools["ceil_div"]

N = 1 << 15
M = 1 << 34
H = 256
L = 6
V = 6**L

def one_digit_tables():
    # cos(2*pi*j/N)
    cos = math_tools["cosines"](N)
    tables = []
    for d in range(3):
        D = set(range(6))-{d}
        counts = [sum((v+h in D) for v in D) for h in range(1,6)]
        table = []
        for i in range(N):
            s = 5*SCALE + 2*sum(c*cos[min((i*h)%N,(-i*h)%N)]
                                for h,c in enumerate(counts,1))
            t = ceil_div(s*M*M,25*SCALE)
            z = isqrt(t)
            table.append(min(M,z+(z*z<t)))
        tables.append(table)
    return tables

def square_upper(table):
    caps = {2:Q(0),3:Q(0)}
    # Reflection reduces 0 <= x <= 1 to 0 <= x <= 1/2.
    for h in range(H//2):
        products = [M]*V
        v = V
        while v>1:
            # f(t) is 12-Lipschitz by centering the digits at 5/2:
            # (2*pi/5) * sum |j-5/2| <= (2*pi/5)*(17/2) < 12.
            guard = ceil_div(6*M,N)+ceil_div(6*M,H*v)
            local = []
            for i in range(v):
                midpoint_num = i*2*H+2*h+1
                near = ((N*midpoint_num+H*v)//(2*H*v))%N
                local.append(min(M,table[near]+guard))
            for i in range(V):
                products[i] = ceil_div(products[i]*local[i%v],M)
            v//=6
        for p in caps:
            J = (6//p)**L
            sums = [0]*J
            for i,f in enumerate(products):
                sums[i%J]+=f
            val = Q(sum(s*s for s in sums),J*M*M)
            caps[p] = max(caps[p],val)
    return caps

def main():
    if not __debug__:
        raise RuntimeError("Assertions are required")
    for digit,table in enumerate(one_digit_tables()):
        caps = square_upper(table)
        for p,tau in ((2,Q("-.300")),(3,Q(".135"))):
            c=caps[p]
            print("d",digit,"p",p,"upper block cap",float(c),flush=True)
            assert c**tau.denominator < Q(6)**(L*tau.numerator)
    print("C_2(L) << 6^(-.300 L), C_3(L) << 6^(.135 L)")

if __name__=="__main__":
    main()
