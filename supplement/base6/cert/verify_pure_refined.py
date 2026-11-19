"""Adaptive rational enclosures of pure-prime grid moments in radix six."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

basic = run_path(str(Path(__file__).with_name("verify_moments.py")))
SCALE = basic['WEIGHT_SCALE']
ceil_div = basic['ceil_div']
TABLE_SIZE = 65536
MOMENTS = (Q(1), Q(4, 3), Q(3, 2), Q(7, 4))
CASES = {
    2: (7, (Q('.108'), Q('.052'), Q('.038'), Q('.025'))),
    3: (6, (Q('.207'), Q('.096'), Q('.067'), Q('.040'))),
}


def rational_cap(modulus, exponent):
    # Largest integer c <= SCALE * modulus**exponent, rounded downward.
    a, b = exponent.numerator, exponent.denominator
    target = SCALE**b * modulus**a
    lo, hi = 0, modulus*SCALE+1
    while hi-lo > 1:
        mid = (lo+hi)//2
        if mid**b <= target:
            lo = mid
        else:
            hi = mid
    return lo


def upper_power_table(square, moment):
    a, b = moment.numerator, 2*moment.denominator
    units = SCALE**(b-a)
    return [basic['upper_root'](x**a*units, b) for x in square]


def certify(prime, degree, power, moment, exponent):
    modulus = prime**degree
    cap = rational_cap(modulus, exponent)
    # Labels h are centers of the cells [(h-1)/n,(h+1)/n].
    # Reflection of the digit amplitudes makes shifts symmetric about 1/2.
    todo = [(1, 4)]
    checks = 0
    max_depth = 0
    frequencies = [6**j for j in range(degree)]
    while todo:
        h, n = todo.pop()
        checks += 1
        assert n <= 1 << 28, ('increase precision or relax the bound', h, n, checks)
        max_depth = max(max_depth,n.bit_length()-1)
        denominator = n*modulus
        errors = [ceil_div(11*moment.numerator*SCALE*(2*TABLE_SIZE*f+denominator),
                           moment.denominator*2*TABLE_SIZE*denominator)
                  for f in frequencies]
        total = 0
        for a in range(modulus):
            phase = ((a*n+h)*TABLE_SIZE) % (TABLE_SIZE*denominator)
            prod = SCALE
            for err in errors:
                k = ((phase+denominator//2)//denominator) % TABLE_SIZE
                prod = ceil_div(prod*min(SCALE, power[k]+err), SCALE)
                phase = (6*phase) % (TABLE_SIZE*denominator)
            total += prod
            if total > cap:
                break
        if total > cap:
            todo.extend([(2*h-1, 2*n), (2*h+1, 2*n)])
    return checks,max_depth


def main():
    if not __debug__:
        raise RuntimeError('Use assertions for all comparisons.')
    cosine = basic['cosines'](TABLE_SIZE)
    for d in (2, 1, 0):
        squared = basic['squared_weights'](d, cosine, TABLE_SIZE)
        for i, moment in enumerate(MOMENTS):
            weight = upper_power_table(squared, moment)
            for prime, (degree, exponents) in CASES.items():
                exponent = exponents[i]
                checks, depth = certify(prime, degree, weight, moment, exponent)
                print(f'd={d}/{5-d}, p={prime}, s={moment}, exponent={exponent}, '
                      f'cells={checks}, final depth={depth}', flush=True)


if __name__ == '__main__':
    main()
