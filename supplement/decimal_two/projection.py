"""Rational tests for the multiplicative projection at exponent .770."""

from fractions import Fraction as Q

from affine import delta, kap2_range, kap5_range, log_bounds


def verify():
    lo2, hi2 = log_bounds(2)
    lo54, hi54 = log_bounds(5, 4)
    lo_k2 = lo2 / (3 * hi2 + hi54)
    hi_k2 = hi2 / (3 * lo2 + lo54)
    assert kap2_range[0] < lo_k2 < hi_k2 < kap2_range[1]
    assert delta > 1 - 3 * lo_k2
    assert Q('.903') < 3 * lo_k2

    # The linear Postnikov fibre gives pointwise decay q^(-.201).
    a, a2, a5 = map(Q, ('.428', '.028', '.2241'))
    pointwise = Q('.201')
    for prime_log, other_log, pure in (
        (kap2_range, kap5_range, a2), (kap5_range, kap2_range, a5)
    ):
        assert (1 - delta - pure) * prime_log[0] / 2 > (
            pointwise * (prime_log[1] + other_log[1] / 2)
        )
    assert (1 - delta - a) / 2 > pointwise
    gamma, y = Q('.565'), Q('.300')
    assert -y + gamma * (Q('2.4') * y - pointwise) < Q('-.0067')
    assert -y + gamma * (Q('1.2') * y + Q('.0485')) < Q('-.069')
    assert Q('1.2') * gamma < 1
    assert Q('2.4') * (Q('.010') + Q('.001')) < 1

    # Cubic and quartic fibre estimates with the directed moment exponents.
    # D_0^delta = q^delta Q^(-delta) supplies the remaining cubic power.
    assert 1 - 2 * Q('.225') == Q('.55')
    # For the quartic power, use 3c + delta = 1 + 2c and c > .903.
    assert 1 + 2 * Q('.903') - 3 * Q('.285') == Q('1.951')
    eta = Q('.770')
    for gamma in (Q('.565'), Q('.665')):
        constant = 1 + Q('1.849') * gamma + Q('.55') * eta / 6
        assert (constant + Q('.006')) ** 2 < 8 * gamma
    for gamma in (Q('.665'), eta):
        constant = 1 + Q('1.7865') * gamma + (Q('.097') / 12 + Q('1.951') / 6) * eta
        assert (constant + Q('.00025')) ** 2 < 9 * gamma
    print('Projection: exact density inequalities through .770 passed.', flush=True)
