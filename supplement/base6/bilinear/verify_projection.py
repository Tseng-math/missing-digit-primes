"""Exact scalar checks for low_projection.md. Run with Python assertions enabled."""
from fractions import Fraction as F

if not __debug__:
    raise RuntimeError("Assertions are required")

eta = F(".76")
delta = F(".102")
a = F(".22947")
b = F(".28965")
assert 5**500 > 6**449
assert 2**500 > 12**139
assert 1+2*(1-delta)-3*b > F("1.925")
assert (delta - (1-2*a)*F(".278"))/3 < -F(".016")

def below(x, z, ceiling):
    """Check (sqrt(x)-1)^2 + z < ceiling by rational squaring."""
    pos = 1+x+z-ceiling
    assert pos > 0 and pos**2 < 4*x

for h in (".5", ".52"):
    g=F(h)
    below(2*g,-F(".016")*g,-F(".0075"))
for h in (".52", ".67"):
    g=F(h)
    below(2*g,(delta*g-(1-2*a)*(g-eta/2))/3,-F(".0045"))
for h in (".67", ".76"):
    g=F(h)
    below(9*g/4,(delta*(g+eta/3)-F("1.925")*(g-2*eta/3))/4,-F(".0014"))
print("Projection exponent 19/25 verified.")
