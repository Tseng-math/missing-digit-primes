"""Exact checks for the orbit estimates in base7_orbits.tex. No dependencies."""
from fractions import Fraction as F
from itertools import combinations


def vertices():
    # A row (v,t) denotes the hyperplane v.x=t for x=(s,b_s,b_c,v).
    unit_rows = [tuple(F(int(i == j)) for i in range(4)) for j in range(4)]
    budget = tuple(map(F, (".249", ".806", ".806", ".903")))
    cuts = [(row, F(0)) for row in unit_rows]
    cuts.extend([(budget, F(".0477")), (unit_rows[0], F(".117"))])
    # .245(e+b_s)-b_s=.92*b_s for each of the two affine pieces of e.
    cuts.extend([((F(0), F("1.675"), F(0), F(0)), F(".245") * F(".117")),
                 ((F("-.245"), F("1.675"), F(0), F(0)), F(0))])
    result = set()
    for rows in combinations(cuts, 4):
        matrix = [list(row) + [rhs] for row, rhs in rows]
        for j in range(4):
            pivot = next((i for i in range(j, 4) if matrix[i][j]), None)
            if pivot is None:
                break
            matrix[pivot], matrix[j] = matrix[j], matrix[pivot]
            matrix[j] = [x / matrix[j][j] for x in matrix[j]]
            for i in range(4):
                if i != j:
                    factor = matrix[i][j]
                    matrix[i] = [a - factor*b for a,b in zip(matrix[i], matrix[j])]
        else:
            x = tuple(matrix[i][-1] for i in range(4))
            if min(x) >= 0 and sum(a*b for a,b in zip(budget,x)) <= F(".0477"):
                result.add(x)
    return result


def values(x):
    s, bs, bc, v = x
    e = max(s, F(".117"))
    g = e + bs
    w = F(".428")*(bs+v) + F(".16")*bc
    ell = F(".0477")-F(".249")*s-F(".806")*(bs+bc)+F(".097")*v
    z = F(".48847")
    tau = z + max(ell,w) + F(".00032")
    y, density = F("0.2"), F(".48")
    rhs = -y + density*(g+bc) + w - bc/2 - max(F(".245")*g-bs,F(".92")*bs)
    return {
        "conductor exponent": g+bc,
        "initial block separation": e+4*(bs+bc)+v,
        "first mass": w,
        "conductor plus height": g+bc+tau,
        "point low curvature": rhs+density*F(".230"),
        "point high curvature": rhs+(density-F(".245"))*tau,
        "square at one fifth": -y+density*(g+bc+tau)/2+F(".04")*(g+v+z)+F(".16")*bc,
        "no-zero square": -F(".5")+(g+bc+tau)/2+F(".04")*(g+v+z)+F(".16")*bc,
    }


def main():
    points = vertices()
    cutoffs = {"conductor exponent": ".192", "initial block separation": ".354",
               "first mass": ".026", "conductor plus height": ".80",
               "point low curvature": "-.003", "point high curvature": "-.003",
               "square at one fifth": "-.003", "no-zero square": "-.06"}
    for key, bound in cutoffs.items():
        best = max(values(x)[key] for x in points)
        if not best < F(bound):
            raise ArithmeticError((key,best,bound))
        print(key, best, "<", bound)


if __name__ == "__main__":
    main()
