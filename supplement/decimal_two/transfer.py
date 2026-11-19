"""Directed positive transfer certificate for normalized decimal word moments."""

from fractions import Fraction as Q

from exact_trig import cos_interval, SCALE

WIDTH = 1 << 26
HEIGHT = 1 << 30


def root(n, m, guess):
    """Floor of n^(1/m); use only verified integer decisions."""
    if guess**m <= n < (guess+1)**m:
        return guess
    lo, hi = 0, WIDTH+1
    while hi-lo > 1:
        mid = (lo+hi)//2
        if mid**m <= n:
            lo = mid
        else:
            hi = mid
    return lo


def weights(digits, cos_nodes):
    M = 2*(len(cos_nodes)-1)
    counts = [sum(d+r in digits for d in digits) for r in range(1,10)]
    out = []
    for t in range(M):
        numerator = 8*SCALE
        for r, count in enumerate(counts, 1):
            j = (r*t)%M
            numerator += 2*count*cos_nodes[min(j,M-j)]
        u = numerator*WIDTH
        q = -(-u//(64*SCALE))
        out.append(min(WIDTH,max(0,q)))
    return out


def certify(digits, cos_nodes, p, lam):
    N = (len(cos_nodes)-1)//5
    assert N%10 == 0
    a,b = p.numerator,p.denominator
    assert 1 <= p <= 2

    def power(w):
        n = w**a * WIDTH**(2*b-a)
        # The exact inequality in root() validates this initial guess.
        guess = int(WIDTH*(w/WIDTH)**(a/(2*b)))
        y = root(n,2*b,guess)
        return y + (y**(2*b) < n)

    powers = [power(w) for w in weights(digits, cos_nodes)]

    def apply(v):
        out = []
        for j in range(N):
            i, rem = divmod(j,10)
            value = 0
            for r in range(10):
                k = r*(N//10)+i
                value += powers[r*N+j]*((10-rem)*v[k]+rem*v[(k+1)%N])
            out.append(value)
        return out

    v = [HEIGHT]*N
    # Find positive approximate eigenvector of the discretized operator.
    # Its accuracy is irrelevant to soundness: check every node below.
    for _ in range(22):
        t=apply(v)
        maximum=max(t)
        v=[HEIGHT*z//maximum for z in t]
    assert min(v)>0 and max(v)<=HEIGHT

    maxdiff = max(abs(v[(j+1)%N]-v[j]) for j in range(N))
    h_lip = Q(N*maxdiff,HEIGHT)
    diff = min(lam.numerator*10*WIDTH*h-lam.denominator*t
               for t,h in zip(apply(v),v))
    nodal = Q(diff,lam.denominator*10*WIDTH*HEIGHT)
    # The distributional second derivative of T_p h is bounded below by
    # -p(80+4 Lip(h)) on each cell. Apply the parabolic chord estimate.
    residual = nodal-p*(80+4*h_lip)/(8*N*N)
    assert residual>0, (digits, p, lam, N, residual)
    return residual
