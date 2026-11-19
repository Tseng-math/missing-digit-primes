"""Rational base-seven endpoint comparisons."""
from fractions import Fraction as F


def check_endpoints():
    assert 7**7921 * 6**100000 < 7**100000 < 7**7923 * 6**100000
    delta_lo, delta_hi = F('.07921'), F('.07923')
    h_min = (1-2*F('.19684'))/2-delta_hi
    h_max = (1-2*F('.19684'))/2-delta_lo
    assert F('.09')*delta_lo > F('.007')
    assert (delta_lo-F('.016'))/3-F('.06')*h_max > F('.007')
    assert h_min+F('.0035') >= F('.22743')
    gamma = F('.770')
    assert (1+2*gamma-gamma*F('.22743')/3+F('.00030'))**2 < 8*gamma
    assert 7**9207 < 6**10000 < 7**9209
    assert (2*F(".9207") - 1 - 2*F(".19684"))/6 > F(".0746")
    assert -F(30, 77) + F(6, 13) - F(".0746") < -F(".0026")
    a = F(".401")
    for m in map(F, (".4195", ".4801")):
        n = 1-m
        lam = m+F(".09696")
        assert m+F(".09693") < lam
        assert F(".08")*(1-lam) < F(".0387")
        for s in (F(".220"), lam):
            exponent = F(".08")*(1-s) - F(".289")*s
            exponent += max(0, s-m) + max(0, s-n)/2
            assert exponent <= -F(".00118")
        for s in (F(".234"), lam):
            correction = F(".0387")+F(".192")*a - (1-a)*s
            correction += max(0,s-n)+max(0,s-m)
            assert correction <= -F(".02447")
    assert F(".0387")-F(".289")*F(".192") < -F(".016")
    assert F(".19389") < F(".4195")-F(".0001")-F(".220")
    y, gamma, z = map(F, (".234", ".192", ".4836"))
    density = F("2.4")*y
    assert -y+gamma*(density-F(".201"))+F(".230")*density < -F(".03559")
    assert -y+gamma*(density-F(".201"))+z*(density-F(".2375")) < -F(".00803")
    assert -y+(F(".0485")+density/2)*(gamma+z) < -F(".01152")
    assert -F(".5")+F(".5485")*(gamma+z) < -F(".1294")
    print("Projection, row-square and small-period endpoint inequalities passed.", flush=True)


if __name__ == '__main__':
    check_endpoints()
