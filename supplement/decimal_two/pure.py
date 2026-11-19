"""Exact uniform shifted first and square norms on pure 2- and 5-grids."""

from fractions import Fraction as Q
from itertools import combinations
from math import isqrt
from subprocess import run

from build import executable
from exact_trig import SCALE, cos_interval
from moments import check_power
from pure_squares import REPRESENTATIVES, squared_norm


NODES = 1 << 17
UNIT = 1 << 36

# Two-adic block length, first cap at 2, first cap at 5 (length eight).
FIRST_PROFILES = {
    (0, 1): (12, '1.008', '7.6'),
    (0, 2): (20, '1.405', '11.7'),
    (0, 3): (12, '1.003', '12.5'),
    (0, 4): (20, '1.43', '13.2'),
    (1, 2): (12, '1.003', '11.2'),
    (1, 3): (20, '1.41', '14.2'),
    (1, 4): (12, '1.008', '14.6'),
    (1, 5): (20, '1.42', '16'),
    (1, 7): (20, '1.42', '14'),
    (1, 8): (12, '1.008', '9.8'),
    (2, 3): (12, '1.003', '11.8'),
    (2, 4): (20, '1.47', '17.9'),
    (2, 5): (12, '1.008', '15.1'),
    (2, 6): (24, '1.585', '14.5'),
    (3, 4): (12, '1.005', '12'),
    (3, 5): (20, '1.45', '16.8'),
    (3, 6): (12, '1.008', '16.8'),
    (4, 5): (12, '1.008', '8.8'),
}


def digit_set(pair):
    return set(range(10)) - set(pair)


def check_alphabets():
    expected = {pair for pair in combinations(range(10), 2) if pair <= (9 - pair[1], 9 - pair[0])}
    covered = set(FIRST_PROFILES)
    assert covered == set(REPRESENTATIVES)

    # The transfer identity in the paper applies to each of these three
    # verified factorizations. It loses a bounded constant for the whole word.
    for left, right in (((0, 5), (0, 1)), ((1, 6), (2, 3)), ((2, 7), (4, 5))):
        base = set(range(5)) - {left[0]}
        assert digit_set(left) == {s + 5 * e for s in base for e in (0, 1)}
        assert digit_set(right) == {2 * s + e for s in base for e in (0, 1)}
        assert right in covered
        covered.add(left)

    for left, right in (((0, 6), (0, 4)), ((0, 7), (0, 3)), ((0, 8), (0, 2))):
        assert digit_set(left) == {10 - d for d in digit_set(right)}
        covered.add(left)
    assert digit_set((0, 9)) == {d - 1 for d in digit_set((0, 1))}
    covered.add((0, 9))
    assert covered == expected


def upper_amplitudes(pair, cosine):
    digits = digit_set(pair)
    counts = [sum(d + h in digits for d in digits) for h in range(1, 10)]
    amplitudes = []
    for j in range(NODES // 2 + 1):
        upper_square = 8 * SCALE
        for h, count in enumerate(counts, 1):
            pos = (h * j) % NODES
            upper_square += 2 * count * cosine[min(pos, NODES - pos)]
        assert upper_square >= 0
        norm_square = (upper_square * UNIT * UNIT + 64 * SCALE - 1) // (64 * SCALE)
        value = isqrt(norm_square)
        amplitudes.append(min(UNIT, value + (value * value < norm_square)))
    return ' '.join(map(str, amplitudes + amplitudes[-2:0:-1]))


def verify_first():
    cosine = [cos_interval(2 * j, NODES)[1] for j in range(NODES // 2 + 1)]
    with executable('pure_first.cpp') as (binary, env):
        for pair, (length_2, cap_2, cap_5) in FIRST_PROFILES.items():
            amplitudes = upper_amplitudes(pair, cosine)
            for prime, length, cap, exponent in (
                (2, length_2, Q(cap_2), Q('.028')),
                (5, 8, Q(cap_5), Q('.2241')),
            ):
                check_power(prime, length * exponent, cap)
                result = run(
                    [str(binary)], input=f'{prime} {length} {cap.numerator} {cap.denominator}\n{amplitudes}\n',
                    text=True, check=True, capture_output=True, cwd=binary.parent, env=env,
                )
                print(f'{pair} {result.stdout.strip()}', flush=True)


def verify_squares():
    for prime, exponent in ((2, Q('.004')), (5, Q('.0053'))):
        for pair in REPRESENTATIVES:
            if prime == 5:
                length, cap = 24, Q('1.223')
            elif pair == (2, 6):
                length, cap = 30, Q('1.084')
            else:
                length, cap = 24, Q('1.068')
            check_power(prime, length * exponent, cap)
            assert squared_norm(*pair, prime, length) < cap, (pair, prime)
        print(f'Pure {prime}-grid: all square witnesses give exponent {exponent}.', flush=True)


def verify():
    check_alphabets()
    verify_squares()
    verify_first()
