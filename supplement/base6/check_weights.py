"""Check the Richert coefficient and parameter reductions for Theorem six:rough."""

from fractions import Fraction as Q


def log_interval(a):
    t = (a - 1) / (a + 1)
    partial = 2 * sum(t**(2*i+1) / (2*i+1) for i in range(30))
    tail = 2*t**61 / (61*(1-t*t))
    return partial, partial+tail


level, rough, mark = Q(14, 25), Q(7, 50), Q(49, 100)
penalty = 1 / (3 - 1 / mark)
assert level / rough == 4 and rough < mark < level < 1 - Q('.43389')
assert (level - rough) / 3 <= rough and 1 - penalty*(1 - rough/mark) > 0
assert 6**43 < 5**50
lower_sieve = log_interval(level/rough - 1)[0]
upper_integral = (log_interval(mark/rough)[1] - (level/mark - 1) *
                  log_interval((level - rough)/(level - mark))[0])
constant_lower = 2/level * (lower_sieve - penalty*upper_integral)
if constant_lower <= Q(212, 1000):
    raise ArithmeticError("The claimed lower bound was not established.")
print("Proved that the two-almost-prime coefficient exceeds 0.212.")
