"""Six-adic orbit exponents. Rational affine vertex enumeration."""
from itertools import product, combinations
from fractions import Fraction as Q

def linear(const=0, *coefficients):
    return (Q(const),) + tuple(map(Q, coefficients or (0, 0, 0, 0)))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def multiply(value, affine):
    return tuple(Q(value) * x for x in affine)


def subtract(a, b):
    return add(a, multiply(-1, b))


def evaluate(affine, point):
    return sum(a * b for a, b in zip(affine, (Q(1), *point)))


def intersection(equations):
    matrix = [list(a[1:]) + [-a[0]] for a in equations]
    for k in range(4):
        pivots = [row for row in range(k, 4) if matrix[row][k]]
        if not pivots:
            return None
        pivot = pivots[0]
        matrix[k], matrix[pivot] = matrix[pivot], matrix[k]
        r = matrix[k][k]
        matrix[k] = [x / r for x in matrix[k]]
        for row in range(4):
            if row != k and matrix[row][k]:
                r = matrix[row][k]
                matrix[row] = [x - r * y for x, y in zip(matrix[row], matrix[k])]
    return tuple(row[-1] for row in matrix)


def vertices(inequalities):
    points = set()
    for equations in combinations(inequalities, 4):
        point = intersection(equations)
        if point is not None and all(evaluate(a, point) >= 0 for a in inequalities):
            points.add(point)
    return points


def admissible_mixture(values, margin):
    """Find a convex mixture bounded by -margin at every vertex."""
    lo, hi = Q(0), Q(1)
    for point, square in values:
        delta = point - square
        rhs = -margin - square
        if delta > 0:
            hi = min(hi, rhs / delta)
        elif delta < 0:
            lo = max(lo, rhs / delta)
        else:
            assert rhs >= 0, (point, margin)
    assert lo <= hi, (lo, hi, margin)


a, delta, sigma, budget, orbit_floor = map(Q,('.434','.102','.2395','.0523','.092'))
z, u_max = Q('.5101'), Q('.240')
s = linear(0,1,0,0,0)
bs = linear(0,0,1,0,0)
bc = linear(0,0,0,1,0)
v = linear(0,0,0,0,1)
weight = linear(0,0,a,2*delta,a)
offset = linear(budget,-sigma,2*delta-1,2*delta-1,delta)
allowed = (s,bs,bc,v,subtract(offset,v))


def suffix(branches):
    eidx,tidx,gidx,kidx=branches
    es=(s,linear(orbit_floor))
    e=es[eidx]
    tt=(offset,weight)
    t=tt[tidx]
    gg=(bs,multiply('.386',e))
    g=gg[gidx]
    gam=add(e,bs)
    kk=(subtract(multiply('.208',gam),bs),add(multiply('.5',bs),multiply('.398',g)))
    saving=kk[kidx]
    constraints=(subtract(e,es[1-eidx]),subtract(t,tt[1-tidx]),subtract(gg[1-gidx],g),subtract(saving,kk[1-kidx]))
    pts=vertices((*allowed,*constraints))
    tt=add(t,linear(z+Q('.00032')))
    d=add(gam,bc)
    qT=add(d,tt)
    pcost=subtract(subtract(weight,multiply('.5',bc)),saving)
    pconst=subtract(pcost,multiply('.232',tt))
    qconst=add(multiply(delta/2,add(add(gam,v),linear(z))),multiply(2*delta,bc))
    return [(p,evaluate(qT,p),evaluate(pconst,p),evaluate(qconst,p),evaluate(d,p),evaluate(pcost,p)) for p in pts]


def run():
    if not __debug__:
        raise RuntimeError('Run without -O')
    data = [v for t in product(range(2), repeat=4) for v in suffix(t)]
    eps = Q('.0003')
    print('Affine vertex records',len(data))
    assert max(v[1] for v in data)<Q('.729')
    assert sigma <= 3*Q('.898')/4 - a
    assert sigma < (1-2*Q('.151'))*Q('.386')
    assert sigma < (1-2*Q('.231'))*Q('.613')
    assert Q('.386')*(Q('.898')-Q('.151'))/(1+Q('.386'))>Q('.208')
    assert Q('.613')*(Q('.898')-Q('.231'))/(1+Q('.613'))>Q('.208')
    assert (Q('.898')-a)/2>=Q('.232')
    y0, y1 = Q('.18'), Q('.24')
    for _,size,hp,sq,cond,cost in data:
        # Before the fixed positive deficit use the log-loss Huxley bound.
        assert -y0 + Q('2.401')*y0*size + hp < -eps
        for y in (y0,y1):
            dens = Q('2.401')*y
            assert -y + dens*(cond+u_max)+max(Q(0),dens-(1-a))*u_max + cost < -eps
            assert -y + dens*(cond+Q('.07')) + cost < -eps
        # For large deficits the high-curvature square exponent decreases.
        assert -y1+Q(15,13)*y1*size+sq < -eps
        # For heights too low for the T-smooth inequality use Huxley's slope.
        assert -y1+Q('1.2005')*y1*(cond+Q('.07'))+sq < -eps
        assert -Q('.5')+size/2+sq < Q('-.09')
    for j in range(1800,2400):
        left, right = Q(j,10000), Q(j+1,10000)
        vals=((-left+Q(30,13)*right*size+hp,
               -left+Q(15,13)*right*size+sq)
              for _,size,hp,sq,cond,cost in data)
        admissible_mixture(vals,eps)
    # Geometric and arithmetic restrictions for mixed denominators.
    a_upper = Q('.102')*Q('.5101')
    assert a_upper +Q('.00021') < budget
    smax = budget/sigma
    bmax = budget/(1-2*delta)
    assert smax < Q('.219')
    assert bmax < Q('.066') and bmax < 2*Q('.386')*orbit_floor
    assert max(smax,orbit_floor+4*bmax) < Q('.356') < Q('.4228')
    assert Q('.4228')-Q('1.586')*Q('.1022')>smax
    assert Q('.4228')-Q('.632')*Q('.1022')>smax
    # Small-order overflow, with log K <= s and s<=.300.
    m0, s0 = Q('.4228'),Q('.305')
    assert a_upper-(Q('.5')-delta)*m0+(Q('.5')-sigma)*s0+Q('.0001')<Q('-.01')
    assert delta-m0+(1-delta-sigma)*s0<0
    assert delta*(1-s0)-sigma*s0<0
    # Partly signed correction at an orbit E smaller than the usual lift.
    wmax=a*bmax
    A = a_upper
    c=Q('.898')
    for m,h in product(map(Q,('.4229','.4801')),map(Q,('.067','.1023'))):
        for e in (orbit_floor,smax):
            for xi in (2*e,m+h):
                if xi>m+h: continue
                corr = (A+wmax-(1-a)*xi+(2-a-3*c/2)*e
                        +max(0,xi-m)+max(0,xi-(1-m))/2)
                assert corr < -Q('.002'), (m,h,e,xi,corr)
    print('600 deficit intervals passed, uniform exponent margin >=',eps)
    print('Small-order, orbit fitting and shortened-lift bounds passed')

if __name__=='__main__':
    run()
