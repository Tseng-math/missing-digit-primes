"""Exact exponent checks for the grouped comparison and its sieve application."""

from fractions import Fraction as F
from itertools import product


DELTA = F(".0643")
ETA_MIN = F(".740")
LAMBDA_MIN = F(".589")
GAMMA_OVERLAP = F(".0002")
RATIONAL_LIMIT = F(".037")
SMALL_SMOOTHING = F(".009")
SMALL_PADDING = F(".012")
SMALL_PERIOD = F(".120")
HIGH_G_FLOOR = F(".34")
HEIGHT_MIN = F(".065")
INITIAL_LOWER_BAND = F(".119")
DIRECT_RESERVE = F(".001")


def check_analytic_range(a, nu_min, height, *, smoothing, padding):
    slope = 1 - 2 * a
    transition = 1 + height - LAMBDA_MIN
    nu_corners = (nu_min, transition, F(1, 2))
    assert nu_min < transition < F(1, 2)
    maximum = {}

    def save(name, value):
        assert value < 0, (name, value)
        maximum[name] = max(value, maximum.get(name, value))

    # On either side of the transition all tests below are affine in
    # nu and eta. A minimum in a denominator threshold is checked by
    # substituting each alternative threshold separately.
    for nu, eta in product(nu_corners, (ETA_MIN, F(1))):
        mu = 1 - nu
        cut = max(LAMBDA_MIN, mu + height)
        largest_row_window = cut - mu
        upper_density = DELTA * (eta - cut)
        gamma_order = F(".002") + upper_density / slope
        small_offset = 1 - cut + SMALL_PADDING
        high_offset = nu + GAMMA_OVERLAP
        assert ETA_MIN > cut > mu
        assert 1 - cut + padding < nu - RATIONAL_LIMIT
        assert 2 * (gamma_order + RATIONAL_LIMIT) < mu - GAMMA_OVERLAP
        assert gamma_order + 3 * RATIONAL_LIMIT + largest_row_window < mu - GAMMA_OVERLAP

        save("diagonal", DELTA * eta - HEIGHT_MIN)
        save("untwisted, capped", -slope * INITIAL_LOWER_BAND + upper_density)
        save("untwisted, uncapped", DELTA * eta - (slope + DELTA) * (nu + HEIGHT_MIN))
        save("nonhorizontal, capped (lower endpoint)", -slope * (mu - GAMMA_OVERLAP) + upper_density)
        save("nonhorizontal, capped (upper endpoint)", 2 * a * cut - mu + upper_density)
        save("nonhorizontal, uncapped", DELTA * eta - (slope + DELTA) * nu + (2 * a - DELTA) * largest_row_window)
        save("horizontal, full order", DELTA * eta - (slope + DELTA) * nu)
        save("horizontal, large order", upper_density - slope * gamma_order)
        save("horizontal, large denominator", upper_density - (1 - 2 * DELTA) * RATIONAL_LIMIT)
        # Dense low-order shells have q/g <= H^2 up to the explicit
        # arbitrarily small completion losses, and g << P.
        assert min(mu - GAMMA_OVERLAP - 2 * HEIGHT_MIN, 2 * mu - cut) > HIGH_G_FLOOR
        save("horizontal, dense", upper_density - slope * HIGH_G_FLOOR + 2 * DELTA * largest_row_window + GAMMA_OVERLAP)

        band_second_moment = (
            DELTA * (gamma_order + largest_row_window)
            + 2 * DELTA * RATIONAL_LIMIT
        )
        save("lifting, boundary smoothing", band_second_moment - smoothing)
        save("lifting, coincident approximant", a * (gamma_order + high_offset) - nu / 2)
        save("lifting, secondary sampling", (a - F(1, 2)) * mu + GAMMA_OVERLAP / 2 + gamma_order / 2)
        assert 1 - 2 * (gamma_order + RATIONAL_LIMIT) - high_offset > 0
        assert 1 - gamma_order - GAMMA_OVERLAP > cut
        assert gamma_order + 3 * RATIONAL_LIMIT + high_offset < 1

        save("small periods, boundary smoothing", DELTA * SMALL_PERIOD - SMALL_SMOOTHING)
        save("small periods, coincident approximant", a * (SMALL_PERIOD + small_offset) - nu / 2)
        save("small periods, secondary sampling", a * mu - F(1, 2) + (SMALL_PERIOD + small_offset) / 2)
        assert small_offset + SMALL_PERIOD < 1

        for secondary in (nu - DIRECT_RESERVE, F(".490"), eta - nu - DIRECT_RESERVE):
            assert secondary > 0
            save("lifting, direct secondary", a * (gamma_order + high_offset + eta - 1) - secondary / 2)
            save("small periods, direct secondary", a * (SMALL_PERIOD + small_offset + eta - 1) - secondary / 2)

    print(f"a={a}, nu_min={nu_min}: least exponent reserve {-max(maximum.values())}")


def verify():
    check_analytic_range(
        F(".3825"), F(".399"), F(".067"), smoothing=F(".020"),
        padding=F(".025"),
    )
    assert DELTA * (F(".1145") + max(F(".089"), RATIONAL_LIMIT)
                    + 2 * RATIONAL_LIMIT) < F(".0180")
    assert F(".0180") < F(".020") < F(".025")
    assert ETA_MIN < F(".741") < F(".742")
    sigma = F(".6175")
    crossover = sigma / (3 - sigma)
    assert 3 * F(".742") / 2 < (1 + crossover)**2 < 3 * F(".742")
    assert -crossover + F(".742") * (sigma + DELTA) / 2 < 0
    bound = F(2, 5)
    width = 1 - 2 * bound
    assert 2 * width >= bound and 2 * (1 - bound) > 1
    assert 1 - bound > F(11, 20) and bound < F(9, 20)
    fee = (3 * bound - 1)**2 / (6 * bound**2 * width)
    assert fee == F(5, 24)
    assert (1 - fee) * F(6, 7) / 7 == F(19, 196) > F(9, 100)
    assert 2 * F(8, 25) == F(".640") < F(".6515")
    word_upper = (F(1, 4) * F(8, 25))**-1
    assert word_upper * F(4, 7) * 7 == 50 < 56
    assert word_upper * F(4, 7) * 7 * F(7, 6) == F(175, 3) < 60
    print(f"Triple loss in the octal proof: {fee}.")


if __name__ == '__main__':
    verify()
