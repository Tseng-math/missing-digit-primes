"""Prove pure-prime block estimates for the mixed radix-6 grids."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

digit_tools = run_path(str(Path(__file__).with_name('verify_moments.py')))
cosines = digit_tools['cosines']
squared_weights = digit_tools['squared_weights']
SCALE = digit_tools['WEIGHT_SCALE']
upper_root = digit_tools['upper_root']

# Rows correspond to the moments 1, 4/3, 3/2 and 7/4.
MOMENTS = (Q(1), Q(4, 3), Q(3, 2), Q(7, 4))
EXPONENTS = {2: (Q('.151'), Q('.081'), Q('.056'), Q('.036')),
             3: (Q('.231'), Q('.119'), Q('.085'), Q('.051'))}
NODES = {2: 3600, 3: 1200}


def four_digit_bounds(prime):
    degree = 4
    modulus = prime**degree
    nodes = NODES[prime]
    size = modulus*nodes
    table = cosines(size)
    for digit in range(3):
        weights = squared_weights(digit, table, size)
        for moment, exponent in zip(MOMENTS, EXPONENTS[prime]):
            a, b = moment.numerator, 2*moment.denominator
            units = SCALE**(b-a)
            powers = [upper_root(x**a*units, b) for x in weights]
            totals = [0]*nodes
            for j, first in enumerate(powers):
                product = first
                index = j
                for _ in range(1, degree):
                    index = (index*6) % size
                    product *= powers[index]
                totals[j % nodes] += product
            guard = moment * Q(137 * sum(36**i for i in range(degree)),
                               8*modulus*nodes**2)
            upper = Q(max(totals), SCALE**degree) + guard
            power = exponent * degree
            assert (upper.numerator**power.denominator <
                    prime**power.numerator * upper.denominator**power.denominator), (prime, digit, moment, float(upper))
            print(f'p={prime}, omitted={digit}/{5-digit}, s={moment}, exponent={exponent}')


def square_bounds(prime):
    degree = 14 if prime == 2 else 11
    modulus = prime**degree
    exponent = Q('.00581') if prime == 2 else Q('.01105')
    occupancy = Q('.0348') if prime == 2 else Q('.0518')
    power = exponent*degree
    for digit in range(3):
        counts = [0]*modulus
        counts[0] = 1
        for j in range(degree):
            shifts = [a * 6**j % modulus for a in range(6) if a != digit]
            counts = [sum(counts[(r-h) % modulus] for h in shifts)
                      for r in range(modulus)]
        upper = Q(modulus*sum(x*x for x in counts), 25**degree)
        assert (upper.numerator**power.denominator <
                prime**power.numerator * upper.denominator**power.denominator)
        print(f'p={prime}, omitted={digit}/{5-digit}, s=2, exponent={exponent}')
        maximum = Q(modulus * max(counts), 5**degree)
        u = degree * occupancy
        assert maximum.numerator**u.denominator < prime**u.numerator * maximum.denominator**u.denominator
        print(f'p={prime}, omitted={digit}/{5-digit}, maximum occupancy exponent={occupancy}')


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError('Use assertions for all comparisons.')
    for p in (2, 3):
        four_digit_bounds(p)
        square_bounds(p)
