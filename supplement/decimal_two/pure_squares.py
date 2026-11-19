"""Exact shifted pure-prime square norm by difference carries."""
from collections import defaultdict
from fractions import Fraction

REPRESENTATIVES = [
    (0,1), (0,2), (0,3), (0,4), (1,2), (1,3), (1,4), (1,5), (1,7),
    (1,8), (2,3), (2,4), (2,5), (2,6), (3,4), (3,5), (3,6), (4,5),
]


def squared_norm(a, b, p, block_length):
    digits = set(range(10)) - {a, b}
    differences = defaultdict(int)
    for d in digits:
        for e in digits:
            differences[d-e] += 1

    modulus = p**block_length
    carries = {0: 1}
    extra_prime = 10//p
    for j in range(block_length):
        next_modulus = modulus//p
        multiplier = pow(extra_prime, j, modulus)
        increments = [[] for _ in range(p)]
        for difference, multiplicity in differences.items():
            step = (multiplier*difference) % modulus
            increments[-step % p].append((step, multiplicity))

        next_carries = defaultdict(int)
        for carry, count in carries.items():
            for step, multiplicity in increments[carry % p]:
                successor = ((carry+step)//p) % next_modulus
                next_carries[successor] += count*multiplicity
        modulus = next_modulus
        carries = next_carries
    assert modulus == 1 and set(carries) == {0}
    return Fraction(p**block_length * carries[0], 64**block_length)
