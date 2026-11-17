#!/usr/bin/env python3
"""Exact rational certificates for all nine base-nine missing-digit factors.

The only input is a fixed positive rational interpolant. All trigonometric
values, fractional powers, residuals, and continuity guards are enclosed
using integer arithmetic. No NumPy/SciPy or floating arithmetic is used.
"""
from fractions import Fraction as F
from math import factorial, isqrt
from pathlib import Path
import hashlib
import json
import time

if not __debug__:
    raise RuntimeError("Run this mathematical certificate without Python -O.")

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "moment_certificates.json"
OUTPUT = ROOT / "moment_verification.json"
S = 10**24


def ceildiv(a, b):
    return -((-a) // b)


def rational(x):
    return f"{x.numerator}/{x.denominator}"


def atan_interval(q, terms):
    a = sum((F((-1)**j, (2*j+1)*q**(2*j+1))
             for j in range(terms)), F(0))
    t = F((-1)**terms, (2*terms+1)*q**(2*terms+1))
    return min(a, a+t), max(a, a+t)


def pi_interval():
    a0, a1 = atan_interval(5, 50)
    b0, b1 = atan_interval(239, 15)
    lo, hi = 16*a0-4*b1, 16*a1-4*b0
    assert F(3) < lo < hi < F(22, 7)
    ilo = lo.numerator*S // lo.denominator
    ihi = ceildiv(hi.numerator*S, hi.denominator)
    assert ihi-ilo <= 1
    return ilo, ihi


def cosine_upper_table(den):
    """Upper bounds S*cos(2*pi*k/den), 0<=k<=den/2."""
    pi_lo, pi_hi = pi_interval()
    terms = 24
    rem = ceildiv(S*4**(2*terms+2), factorial(2*terms+2))
    table = []
    for k in range(den//2+1):
        xlo = pi_lo*(2*k)//den
        xhi = ceildiv(pi_hi*(2*k), den)
        ylo = xlo*xlo//S
        yhi = ceildiv(xhi*xhi, S)
        alo = ahi = slo = shi = S
        for j in range(1, terms+1):
            div = S*(2*j-1)*(2*j)
            alo = alo*ylo//div
            ahi = ceildiv(ahi*yhi, div)
            if j % 2:
                slo -= ahi
                shi -= alo
            else:
                slo += alo
                shi += ahi
        assert slo-rem <= shi+rem
        table.append(min(S, shi+rem))
    return table


def power_upper(v, p, q):
    """Upper bound S*(v/S)^(p/(2q)), for the squared digit modulus v/S."""
    assert 0 <= v <= S
    if (p, q) == (1, 1):
        target = v*S
        r = isqrt(target)
        y = r if r*r == target else r+1
        assert y*y >= target and (y == 0 or (y-1)**2 < target)
    elif (p, q) == (3, 2):
        target = v**3*S
        r = isqrt(isqrt(target))
        y = r if r**4 == target else r+1
        assert y**4 >= target and (y == 0 or (y-1)**4 < target)
    else:
        raise ValueError("Unsupported moment")
    assert 0 <= y <= S
    return y


def verify_one(row, coshi):
    assert row["base"] == 9
    d = row["forbidden_digit"]
    assert d in range(5)
    p, q = row["moment_numerator"], row["moment_denominator"]
    assert (p, q) in ((1, 1), (3, 2))
    s = F(p, q)
    n0, scale0 = row["coarse_grid_size"], row["denominator"]
    refinement, coarse = row["refinement"], row["numerators"]
    assert isinstance(n0, int) and n0 >= 2
    assert isinstance(scale0, int) and scale0 > 0
    assert isinstance(refinement, int) and refinement > 0
    assert len(coarse) == n0
    assert all(isinstance(x, int) and x > 0 for x in coarse)
    g = [(refinement-r)*coarse[i]+r*coarse[(i+1) % n0]
         for i in range(n0) for r in range(refinement)]
    n, gscale = n0*refinement, scale0*refinement
    den = 9*n
    assert len(coshi) == den//2+1
    lam = F(2477, 10000) if s == 1 else F(314229, 2000000)
    counts = [9-r-int(d+r <= 8)-int(d-r >= 0) for r in range(1, 9)]
    assert all(x >= 0 for x in counts)
    worst_num = None
    worst_i = None
    # Common denominator for all residuals avoids Fraction work inside nodes.
    residual_den = 81*S*gscale*lam.denominator
    for i in range(n):
        total = 0
        for j in range(9):
            ell = i+j*n
            square_num = 8*S
            for r, multiplicity in enumerate(counts, 1):
                k = r*ell % den
                k = min(k, den-k)
                square_num += 2*multiplicity*coshi[k]
            square_hi = min(S, max(0, ceildiv(square_num, 64)))
            w = power_upper(square_hi, p, q)
            idx, rem = divmod(ell, 9)
            interpolated = (9-rem)*g[idx]+rem*g[(idx+1) % n]
            total += w*interpolated
        residual_num = total*lam.denominator-81*S*lam.numerator*g[i]
        if worst_num is None or residual_num > worst_num:
            worst_num, worst_i = residual_num, i
    max_residual = F(worst_num, residual_den)
    minimum, maximum = F(min(g), gscale), F(max(g), gscale)
    lip_g = F(n*max(abs(g[(i+1) % n]-g[i]) for i in range(n)), gscale)
    # Center all digit phases at 4; sum_{a=0}^8 |a-4|=20.
    # Lip(f_d) <= (2*pi/8)*20 = 5*pi < 110/7, uniformly in d.
    lip_f = F(110, 7)
    lip_residual = (s*lip_f*maximum+lip_g)/9+lam*lip_g
    guard = lip_residual/(2*n)
    margin = -max_residual-guard
    assert margin > 0, f"Continuum certificate failed for d={d}, s={s}: {margin}"
    return dict(forbidden_digit=d, reflected_digit=8-d,
                moment=rational(s), lambda_bound=rational(lam),
                coarse_grid_size=n0, refinement=refinement,
                verification_grid_size=n, minimum_g=rational(minimum),
                maximum_g=rational(maximum), lipschitz_g=rational(lip_g),
                maximum_node_residual_upper=rational(max_residual),
                maximum_residual_node=worst_i,
                continuous_guard=rational(guard), strict_margin=rational(margin))


def main():
    start = time.monotonic()
    data = json.loads(INPUT.read_text())
    rows = data["rows"]
    assert len(rows) == 10
    assert {(r["forbidden_digit"],r["moment_numerator"],r["moment_denominator"])
            for r in rows} == {(d,p,q) for d in range(5) for p,q in ((1,1),(3,2))}
    assert 2477**200*9**127 < 10000**200
    assert 314229**40000*9**33693 < 2000000**40000
    dimensions = sorted({9*r["coarse_grid_size"]*r["refinement"] for r in rows})
    results = []
    for den in dimensions:
        print(f"Exact cosine enclosure: denominator {den}", flush=True)
        coshi = cosine_upper_table(den)
        for row in rows:
            if 9*row["coarse_grid_size"]*row["refinement"] != den:
                continue
            result = verify_one(row, coshi)
            results.append(result)
            print(f"PASS d={result['forbidden_digit']} s={result['moment']} n={result['verification_grid_size']}", flush=True)
        del coshi
    result = dict(status="PASS: exact integer/rational continuous supersolution certificates",
                  base=9, forbidden_digits=list(range(9)),
                  moment_one_discrete_exponent="73/200",
                  moment_one_continuum_decay="127/200",
                  moment_three_halves_discrete_exponent="6307/40000",
                  moment_three_halves_continuum_decay="33693/40000",
                  kappa="6307/20000",
                  interpolated_four_thirds_discrete_exponent="13607/60000",
                  arithmetic_scale=S, cosine_taylor_degree=48,
                  pi_scaled_interval=pi_interval(),
                  certificate_sha256=hashlib.sha256(INPUT.read_bytes()).hexdigest(),
                  verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  results=results,
                  elapsed_seconds=time.monotonic()-start,
                  scope="Uniform shifted-grid Fourier moments for all nine padded word alphabets; no prime conclusion is asserted by this verifier.")
    # Interpolation: g(4/3) <= (g(1)+2g(3/2))/3.
    assert (F(73,200)+2*F(6307,40000))/3 == F(13607,60000)
    OUTPUT.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k != 'results'}, indent=2))


if __name__ == "__main__":
    main()
