"""Directed interval proof for high difference orders in radix six.

All lengths and exponents, including intermediate enclosures, are in
integer billionths. This program checks a bound for the unmasked Gram.
No input files, dependencies or machine floating-point operations.
"""

U = 10**9
R12 = (12, 16, 18, 21, 24)
MOM = (433890000, 289650000, 229470000, 155800000, 101756000)
PURE = ((108000000, 52000000, 38000000, 25000000, 8130000),
        (207000000, 96000000, 67000000, 40000000, 12150000))
DEL, DELLO = 101756000, 101755000
A = MOM[0]
C = 3*(U-DEL)//4 - A
KLO = (386852000, 613147000)
KHI = (386853000, 613148000)
TAIL = (383000000, 118000000)
M_MIN, M_MAX = 422900000, 480100000
HMIN, HMAX = 77000000, 101856000
TARGET = -100000

def up(x, y):
    return -((-x*y)//U)

def down(x, y):
    return x*y//U

def upper(a, b, d, p, index):
    """a=m, b=s, d=min valuation, p=t0, index names the heavy prime."""
    ml, mr = a
    sl, sr = b
    dl, dr = d
    tl, tr = p
    vl = (ml+102010000)//2
    vr = (mr+102010000+1)//2
    wl, wr = vl-43000000, vr-43000000
    # The top layer is coded by tl=tr=0. The last by tl=tr=1.
    if tl == tr == 0:
        hl = hr = HMAX
        tf = 220000000
    elif tl == tr == U:
        hl = hr = HMIN
        tf = U
    else:
        # t (the actual left shell edge) belongs to [tl,tr].
        # The right edge is at most t+.0005.
        # All rational moment evaluations below are outward for the chosen h.
        def h_bounds(t, lower):
            g = min(v*12+(r-24)*t for v, r in zip(MOM,R12))
            g -= 12*DEL
            num = g*(wr if (g<0)==lower else wl)
            correction = num//(12*U) if lower else -((-num)//(12*U))
            return max(HMIN, DEL+100000+correction)
        hl, hr = h_bounds(tr, True), h_bounds(tl, False)
        tf = tr+500000
    sr = min(sr, mr+hr+30000)
    sl = max(sl,dl)
    dr = min(dr,sr)
    if sl>sr or dl>dr:
        return -U
    k, K, bp = KLO[index], KHI[index], PURE[index]
    # xi=d+k(s-d), increasing in both s and d on the domain d<=s.
    xl = dl+down(k, sl-dl)
    xr = dr+up(K, sr-dr)
    lp = down(k, U-2*bp[0])
    # -C d - L(s-d) = -L s + (L-C)d, with L>=C.
    assert lp >= C
    base = up(DEL, U-sl)-down(lp, sl)+up(lp-C, dr)
    base += max(0, sr-ml)+(max(0, xr-(U-mr))+1)//2
    if xl- mr >= 10000 and hr+xr-ml < dl//2-10000:
        ell = max(0,xl-(U-ml))
        joint = min(down(DELLO,ell)
                    +down(DELLO,min(dl//2,2*(hl+xl-mr)))//2,
                    xl-mr+down(DELLO,ell)//2)
    else:
        joint = 0
    # Lengths of the common filter block lying in the mixed and pure
    # pieces. The v,w endpoint uncertainty is handled monotonically.
    l1 = max(0, min(wl,dl-vr))
    l2 = max(0, min(vl+wl,sl)-max(vr,dr))
    gain = 0
    if p != (U,U):
        for r,f,bpi in zip(R12,MOM,bp):
            penalty = -(-((r-12)*tf*wr)//(12*U))
            gain=max(gain,down(A-f,l1)+down(down(bp[0]-bpi,k),l2)-penalty)
    # Last-L folded-average bound. It is used independently of the
    # joint ratio mean and of the taper gain.
    length = min(max(0,sl-dr), max(0,sl-mr)*U//(U-k))
    tail = down(TAIL[index], length)
    return base-max(joint+2*gain,tail)

def verify(index, pdomain):
    initial = [(M_MIN,M_MAX),(305000000,M_MAX+HMAX+30000),
               (0,M_MAX+HMAX+30000),pdomain]
    stack=[initial]
    checked=0
    cap=-U
    while stack:
        box=stack.pop()
        value=upper(*box,index)
        if value<=TARGET:
            checked+=1
            cap=max(cap,value)
            continue
        # Prefer splitting the three lengths, since their derivatives
        # in the estimate dominate that of the amplitude threshold.
        widths=[(hi-lo)*(1 if i==3 else 4) for i,(lo,hi) in enumerate(box)]
        var=max(range(4),key=widths.__getitem__)
        lo,hi=box[var]
        if hi-lo<=1:
            raise ArithmeticError((index,box,value))
        mid=(lo+hi)//2
        for interval in ((lo,mid),(mid,hi)):
            child=list(box)
            child[var]=interval
            stack.append(child)
    return checked,cap

if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError("Assertions are required.")
    for i in range(2):
        for p in ((0,0),(220000000,500000000),(U,U)):
            bounds=verify(i,p)
            print("heavy",i+2,"shell",p,"boxes",bounds[0],
                  "upper bound in billionths",bounds[1],flush=True)
