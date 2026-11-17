#!/usr/bin/env python3
"""Reconstruct and verify all 188 rational Type II systems.

The analytic constants, windows, margins, variables, branches and objectives
are fixed below. Each stored rational constraint is compared with its exact
reconstruction before the nonnegative dual or Farkas combination is checked.
There are 36 forward systems and 76 systems for each reverse input pair.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from types import SimpleNamespace
import json
import sys

import forward_systems as h
from verify_rational_duals import verify_one

HERE = Path(__file__).resolve().parent
BASE_NAMES = [
    'nu', 'u', 'a0', 'ell', 's', 'j', 'a', 'b', 'c', 'd', 'e', 'f',
    's1', 's2', 'total',
]
COMMON = {
    'generic_margin': F(1, 1000000),
    'arc_margin': F(1, 1000000),
    'width': F(1, 1000000),
}
BASE_PARAMETERS = {
    'nu_min': F(257, 1000),
    'nu_max': F(29, 100),
    **COMMON,
}
STRICT_PARAMETERS = {
    'nu_min': F(249999, 1000000),
    'nu_max': F(129, 500),
    **COMMON,
    'deep_margin': F(1, 1000000),
}
CLASSICAL_P = F(235, 154)
CLASSICAL_G = F(59, 433)
FIXED_GAMMA = F(27, 77)
FIXED_EPS = F(23, 500)
FIXED_G4 = F(3, 5000)
FIXED_KAPPA = F(23, 80)
STRICT_GROUPS = [
    ('central', F(14387, 100000), F(53, 10000)),
    ('uniform', F(287, 2000), F(1, 250)),
]


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def read_record(filename):
    return json.loads((HERE / filename).read_text())


def compare_parameters(actual, expected, label):
    require(set(actual) == set(expected), label + ': parameter keys differ')
    for key, value in expected.items():
        require(F(actual[key]) == value, label + ': wrong ' + key)


def canonical_system(system):
    """Retain labels/order, and compare each numeric entry as a Fraction."""
    require(set(system) in (
        {'names', 'constraints'}, {'names', 'constraints', 'objective'}),
        'unexpected system fields')
    names = tuple(system['names'])
    require(len(set(names)) == len(names), 'duplicate variable names')
    constraints = []
    for row in system['constraints']:
        require(set(row) == {'label', 'coefficients', 'rhs'},
                'unexpected constraint fields')
        require(len(row['coefficients']) == len(names),
                'constraint dimension mismatch')
        constraints.append((
            row['label'], tuple(F(x) for x in row['coefficients']), F(row['rhs'])))
    objective = None
    if 'objective' in system:
        obj = system['objective']
        require(set(obj) == {'coefficients', 'constant'},
                'unexpected objective fields')
        require(len(obj['coefficients']) == len(names),
                'objective dimension mismatch')
        objective = (tuple(F(x) for x in obj['coefficients']), F(obj['constant']))
    return names, tuple(constraints), objective


def verify_rebuilt(item, rows, labels, objective, label):
    require(type(item['feasible']) is bool, label + ': feasibility flag type')
    # An infeasibility certificate has no objective; it must refute the SAME
    # complete constraints, so no optimization or presumption is needed.
    expected = h.encoded_system(
        rows, labels, objective if item['feasible'] else None)
    require(canonical_system(item['system']) == canonical_system(expected),
            label + ': saved system differs from fixed exact reconstruction')
    entries = item['dual'] if item['feasible'] else item['farkas']
    require(len(entries) > 0, label + ': empty certificate')
    for entry in entries:
        require(type(entry['row']) is int and 0 <= entry['row'] < len(rows),
                label + ': invalid dual row index')
    value = verify_one(item)
    require(value < 0, label + ': bound not strictly negative')
    return value


def branch_map(items, bits, label):
    expected = set(product((0, 1), repeat=bits))
    require(len(items) == len(expected), label + ': incorrect branch count')
    result = {}
    for item in items:
        key = tuple(item['branch'])
        require(len(key) == bits and all(type(x) is int and x in (0, 1) for x in key),
                label + ': invalid branch tuple')
        require(key not in result, label + ': duplicated branch')
        result[key] = item
    require(set(result) == expected, label + ': incomplete branch set')
    return result


def verify_base():
    require(h.NAMES == BASE_NAMES and h.D == len(BASE_NAMES),
            'base variable list changed')
    require((h.P, h.G, h.GAMMA, h.EPS, h.G4, h.KAPPA) ==
            (CLASSICAL_P, CLASSICAL_G, FIXED_GAMMA, FIXED_EPS, FIXED_G4, FIXED_KAPPA),
            'base analytic constants changed')
    data = read_record('forward_certificates.json')
    compare_parameters(data['parameters'], BASE_PARAMETERS, 'base')
    expected_basic = {
        'line_trace', 'line_mean_optimized', 'complete_mean',
        'global_complete_mean',
    }
    require(set(data['basic']) == expected_basic, 'base basic objectives differ')
    rows, labels, m = h.setup(**BASE_PARAMETERS)
    ell, s, x, nu, theta = (m[name] for name in ('ell', 's', 'x', 'nu', 'theta'))
    objectives = {
        'line_trace': h.add(ell, h.scale(s, F(3, 2)),
                            h.scale(x, -F(1, 2)), h.scale(nu, -1)),
        'line_mean_optimized': h.add(ell, h.scale(s, 2), h.scale(nu, -1),
                                     h.scale(theta, -F(1, 3))),
        'complete_mean': h.add(ell, nu, h.scale(s, 2), h.scale(theta, -1)),
    }
    values = {}
    for name, objective in objectives.items():
        values[name] = verify_rebuilt(
            data['basic'][name], rows, labels, objective, 'base/' + name)
    global_parameters = dict(BASE_PARAMETERS)
    global_parameters['generic_margin'] = None
    rr, ll, mm = h.setup(**global_parameters)
    objective = h.add(mm['ell'], mm['nu'], h.scale(mm['s'], 2),
                      h.scale(mm['theta'], -1))
    values['global_complete_mean'] = verify_rebuilt(
        data['basic']['global_complete_mean'], rr, ll, objective,
        'base/global_complete_mean')
    branches = branch_map(data['rank2'], 5, 'base')
    rank_values = []
    for branch in sorted(branches):
        rr, ll, mm = h.setup(**BASE_PARAMETERS, branch=branch)
        value = verify_rebuilt(branches[branch], rr, ll, mm['objective'],
                               'base/rank2/' + str(branch))
        if branches[branch]['feasible']:
            rank_values.append(value)
    print('BASE: rebuilt all 36 exact systems; all dual/Farkas checks passed.')
    print('  largest rank2 bound =', max(rank_values))
    print('  line trace bound    =', values['line_trace'])
    return 36


def verify_strict_groups():
    # The reverse system extends the forward variable list by two entries.
    # Import it after verifying the forward systems.
    import reverse_systems as qa
    require(h.NAMES == BASE_NAMES + ['ell0', 'deep'] and
            h.D == len(BASE_NAMES) + 2, 'affine variable list changed')
    require((qa.CLASSICAL_P, qa.CLASSICAL_G) == (CLASSICAL_P, CLASSICAL_G),
            'classical moment parameters differ')
    total_count = 0
    for group, g, delta in STRICT_GROUPS:
        h.P = F(3, 2)
        h.G = g
        h.GAMMA = FIXED_GAMMA
        h.EPS = FIXED_EPS
        h.G4 = FIXED_G4
        h.KAPPA = FIXED_KAPPA
        arguments = {key: str(value) for key, value in STRICT_PARAMETERS.items()}
        arguments.update({
            'p': '3/2', 'g': str(g), 'delta': str(delta),
            'eps': str(FIXED_EPS), 'g4': str(FIXED_G4),
        })
        args = SimpleNamespace(**arguments)
        data = read_record('reverse_' + group + '_certificates.json')
        compare_parameters(data['parameters'], STRICT_PARAMETERS, group)
        compare_parameters(data['analytic_inputs'], {
            'p': F(3, 2), 'g': g, 'delta': delta,
            'eps': FIXED_EPS, 'g4': FIXED_G4,
        }, group + '/analytic inputs')
        expected_basic = {
            name + str(db) for db in (0, 1) for name in (
                'line_trace', 'line_mean', 'base_mean',
                'deep_below_B', 'deep_below_R', 'bad_nonzero')
        }
        require(set(data['basic']) == expected_basic,
                group + ': basic objective set changed')
        values = {}
        for db in (0, 1):
            rows, labels, m = qa.build(args, deep_branch=db)
            objectives = {
                'line_trace': h.add(m['ell'], h.scale(m['s'], F(3, 2)),
                                    h.scale(m['x'], -F(1, 2)), h.scale(m['nu'], -1)),
                'line_mean': h.add(m['ell'], h.scale(m['s'], 2),
                                   h.scale(m['nu'], -1), h.scale(m['theta'], -F(1, 3))),
            }
            for name, objective in objectives.items():
                key = name + str(db)
                values[key] = verify_rebuilt(
                    data['basic'][key], rows, labels, objective, group + '/' + key)
            rr, ll, mm = qa.build(args, deep_branch=db, generic=False)
            objectives = {
                'base_mean': h.add(mm['ell0'], mm['nu'], h.scale(mm['s'], 2),
                                   h.scale(mm['theta'], -1)),
                'deep_below_B': h.sub(mm['deep'], mm['blo']),
                'deep_below_R': h.sub(mm['deep'], mm['nu']),
                'bad_nonzero': h.add(mm['ell0'], mm['s'], h.scale(mm['x'], -1),
                                    h.scale(mm['deep'], -1)),
            }
            for name, objective in objectives.items():
                key = name + str(db)
                values[key] = verify_rebuilt(
                    data['basic'][key], rr, ll, objective, group + '/' + key)
        branches = branch_map(data['rank2'], 6, group)
        rank_values = []
        for branch in sorted(branches):
            db, geometric = branch[0], branch[1:]
            rr, ll, mm = qa.build(args, branch=geometric, deep_branch=db)
            value = verify_rebuilt(branches[branch], rr, ll, mm['objective'],
                                   group + '/rank2/' + str(branch))
            if branches[branch]['feasible']:
                rank_values.append(value)
        print(group.upper() + ': rebuilt all 76 exact systems; '
              'all dual/Farkas checks passed.')
        print('  largest rank2 bound =', max(rank_values))
        print('  largest line trace  =',
              max(values['line_trace0'], values['line_trace1']))
        total_count += 76
    return total_count


def main():
    require(sys.flags.optimize == 0, 'run without -O: exact certificate assertions must remain enabled')
    # Independently check the rational replacements for log_10 9 and g4.
    require(10**477 < 9**500, 'alpha rational bound failed')
    require(730**5000 < 10**3 * 729**5000, 'fourth-moment rational bound failed')
    count = verify_base() + verify_strict_groups()
    require(count == 188, 'wrong total system count')
    prohibited = ('numpy', 'scipy', 'sympy')
    require(not any(name == p or name.startswith(p + '.')
                    for name in sys.modules for p in prohibited),
            'verification imported an optimization/numerical package')
    print('VERIFIED: 188 / 188 freshly rebuilt exact systems.')
    print('All matrix identities and strict rational inequalities verified.')


if __name__ == '__main__':
    main()

