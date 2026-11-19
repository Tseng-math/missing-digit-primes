"""Exact affine-cell checks for long_row_mask_mellin.md."""

from fractions import Fraction as Q
from itertools import combinations, product


def linear(const=0, *coefficients):
    return (Q(const),) + tuple(map(Q, coefficients or (0, 0, 0, 0)))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def multiply(value, affine):
    return tuple(Q(value) * x for x in affine)


def subtract(a, b):
    return add(a, multiply(-1, b))


def evaluate(affine, point):
    return sum(a * b for a, b in zip(affine, (Q(1), *point)))


def intersection(equations):
    matrix = [list(a[1:]) + [-a[0]] for a in equations]
    for k in range(4):
        pivots = [row for row in range(k, 4) if matrix[row][k]]
        if not pivots:
            return None
        pivot = pivots[0]
        matrix[k], matrix[pivot] = matrix[pivot], matrix[k]
        r = matrix[k][k]
        matrix[k] = [x / r for x in matrix[k]]
        for row in range(4):
            if row != k and matrix[row][k]:
                r = matrix[row][k]
                matrix[row] = [x - r * y for x, y in zip(matrix[row], matrix[k])]
    return tuple(row[-1] for row in matrix)


def vertices(inequalities):
    points = set()
    for equations in combinations(inequalities, 4):
        point = intersection(equations)
        if point is not None and all(evaluate(a, point) >= 0 for a in inequalities):
            points.add(point)
    return points


def admissible_mixture(values, margin):
    """Find a convex mixture bounded by -margin at every vertex."""
    lo, hi = Q(0), Q(1)
    for point, square in values:
        delta = point - square
        rhs = -margin - square
        if delta > 0:
            hi = min(hi, rhs / delta)
        elif delta < 0:
            lo = max(lo, rhs / delta)
        else:
            assert rhs >= 0, (point, margin)
    assert lo <= hi, (lo, hi, margin)


def check_cells(*, a, budget, min_short, order_slope, upper_decay, min_orbit,
                linear_saving, cut_at_eta_one, low_break, high_break, height_ceiling):
    floor_height = 1 - min_short - cut_at_eta_one
    tail_reserve = Q('.00032')
    s = linear(0, 1, 0, 0, 0)
    b_s = linear(0, 0, 1, 0, 0)
    b_c = linear(0, 0, 0, 1, 0)
    v = linear(0, 0, 0, 0, 1)
    weight = linear(0, 0, a, '.194', a)
    offset = linear(budget, -order_slope, '-.806', '-.806', '.097')
    allowed = (s, b_s, b_c, v, subtract(offset, v))
    data = []
    for orbit_index, height_index, valuation_index, saving_index in product(range(2), repeat=4):
        orbit_choices = (linear(min_orbit), s)
        orbit = orbit_choices[orbit_index]
        height_choices = (offset, weight)
        height = add(linear(floor_height + tail_reserve), height_choices[height_index])
        valuations = (b_s, multiply('.301', orbit))
        valuation = valuations[valuation_index]
        gamma_s = add(orbit, b_s)
        savings = (subtract(multiply(linear_saving, gamma_s), b_s),
                   add(multiply('.5', b_s), multiply('.403', valuation)))
        saving = savings[saving_index]
        extra = (
            subtract(orbit, orbit_choices[1 - orbit_index]),
            subtract(height_choices[height_index], height_choices[1 - height_index]),
            subtract(valuations[1 - valuation_index], valuation),
            subtract(saving, savings[1 - saving_index]),
        )
        points = vertices((*allowed, *extra))
        if not points:
            continue

        q_t = add(add(gamma_s, b_c), height)
        point_cost = subtract(subtract(weight, multiply('.5', b_c)), saving)
        high_constant = subtract(point_cost, multiply(upper_decay, height))
        lower_size = add(add(gamma_s, b_c), linear('.230'))
        square_constant = add(multiply('.0485', add(add(add(orbit, b_s), v), linear(floor_height))),
                              multiply('.194', b_c))
        for point in points:
            xi = evaluate(q_t, point)
            assert xi < height_ceiling
            high = evaluate(high_constant, point)
            low = evaluate(point_cost, point)
            low_size = evaluate(lower_size, point)
            square = evaluate(square_constant, point)
            data.append((xi, high, low, low_size, square))

    def point_at(y, xi, constant):
        return -y + Q('2.4') * y * xi + constant

    assert data
    assert Q('2.4') * high_break < 1 - a
    assert Q('2.4') * low_break > upper_decay
    assert Q('1.2') * (height_ceiling + Q('.097') * Q('.230')) < 1
    assert Q('1.5') * (height_ceiling + Q('.097') * Q('.230')) / Q('1.25')**2 < 1
    for xi, hi, lo, low_size, sq in data:
        assert xi > 1 / Q('2.4')
        assert xi < 1 / Q('1.2')
        assert point_at(low_break, xi, hi) < Q('-.001')
        assert point_at(low_break, low_size, lo) < Q('-.001')
        assert point_at(high_break, low_size, lo) < Q('-.001')
        assert -high_break + Q('1.2') * high_break * xi + sq < Q('-.0005')

    assert (10000 * low_break).denominator == (10000 * high_break).denominator == 1
    for k in range(int(10000 * low_break), int(10000 * high_break)):
        left, right = Q(k, 10000), Q(k + 1, 10000)
        admissible_mixture(
            ((point_at(right, xi, hi), -left + Q('1.2') * left * xi + sq)
             for xi, hi, _, _, sq in data),
            Q('.0004'),
        )

    return len(data)
