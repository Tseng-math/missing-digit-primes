"""Exact inequalities in the small-period comparison. Decimal strings are rational."""

from fractions import Fraction as F


def main():
    assert 5**1000 > 6**898
    assert 5**1000 < 6**899
    assert 2**46 > 3**27
    assert (F(".898") - F(".43389")) * F(".27") > F(".125")
    assert (F(".898") - F(".43389")) / 2 > F(".232")
    assert 1 - F(".43389") > F(".566")
    g, z, u = F(".235"), F(".479"), F(".241")
    near, cut = F(".04"), F(".215")
    # The first branch maximizes the tau part at tau=u. For y=.04
    # the conductor term and second-branch height term are nonpositive.
    slope = F("2.4") * near
    assert slope < min(F(".125"), F(".232"), F(".566"))
    assert -near + slope * u < F("-.016")
    slope = F(30, 13) * cut
    assert F(".232") < slope < F(".566")
    p_low = -cut + (slope - F(".125")) * g + slope * u
    p_high = -cut + (slope - F(".125")) * g + (slope - F(".232")) * z
    q = -cut + (slope / 2 + F(".051")) * (g + z)
    no_zero = -F(1, 2) + F(".551") * (g + z)
    assert p_low < F("-.008")
    assert p_high < F("-.0012")
    assert q < F("-.0014")
    assert F(15, 13) * (g + z) < 1
    assert no_zero < F("-.1")
    print("All rational inequalities hold.")
    print("Medium-layer upper bounds (exact):", p_low, p_high, q)


if __name__ == "__main__":
    main()
