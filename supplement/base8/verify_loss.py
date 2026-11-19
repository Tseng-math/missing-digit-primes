"""Exact rational enclosure of the limiting octal positive-sieve loss."""

from fractions import Fraction
from math import factorial


def enclosure(degree):
    lower = sum(
        Fraction(factorial(index) * factorial(total - index),
                 2 ** (total + 2) * factorial(total + 2))
        for total in range(degree + 1) for index in range(total + 1)
    )
    return lower, lower + Fraction(1, 2 ** (degree + 3))


if __name__ == '__main__':
    lower, upper = enclosure(30)
    if not (Fraction('.18601165999') < lower <= upper < Fraction('.18601166012')
            and Fraction(6, 49) * (1 - upper) > Fraction('.09967')
            and Fraction(6, 49) * (1 - Fraction(5, 24)) > Fraction(9, 100)):
        raise ArithmeticError('Octal rational enclosure failed')
    print('0.18601165999 < loss < 0.18601166012; real-endpoint liminf > 0.09967')
