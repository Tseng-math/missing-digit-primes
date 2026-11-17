#!/usr/bin/env python3
"""Verify all base-ten moments, Hellinger profiles, and Type II systems.

Run: python3 -B -S scripts/base10/verify_base10.py
All certificate inequalities are checked with exact integers and fractions.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import sys
import time

import verify_three_halves as pressure
import verify_classical_moments as classical_pressure
import verify_hellinger as hellinger
import verify_type_ii as budgets

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    require(sys.flags.optimize == 0, 'Do not use -O: certificate assertions must be enabled.')
    started = time.monotonic()
    require(9**500 > 10**477, 'Lower bound for the digital exponent failed')
    require(9**100 < 10**99, 'Upper bound for the digital exponent failed')
    require(9**21 > 10**20, 'Mixed-grid exponent bound failed')
    require(730**5000 < 10**3 * 729**5000, 'Fourth-moment bound failed')
    require(F(23850, 37879) > F(629, 1000), 'Maximal Type I margin failed')
    print('PART 1: continuous-state classical moment upper certificates', flush=True)
    expected_classical_specs = {
        'one': (F(1), F(27, 77), F(1121, 500), 37, 4),
        'classical': (F(235, 154), F(59, 433), F(1368543, 1000000), 56, 6),
    }
    require(classical_pressure.SPECS == expected_classical_specs,
            'Classical moment parameters changed')
    require(classical_pressure.NODES == 1000 and classical_pressure.DEN == 1 << 40,
            'Classical moment grid or denominator changed')
    classical_results = classical_pressure.verify_all(HERE)
    require([(row['kind'], row['digit']) for row in classical_results] ==
            [(kind, digit) for kind in expected_classical_specs for digit in range(5)],
            'Missing or repeated classical moment certificate')
    print('PART 2: continuous-state 3/2 moment upper certificates', flush=True)
    for digit in range(5):
        name = f'moment_three_halves_d{digit}.json'
        record = json.loads((HERE / name).read_text())
        expected_g = F(14387, 100000) if digit == 4 else F(287, 2000)
        expected_n = 2000 if digit == 4 else 1000
        expected_lambda = F(6963699, 5000000) if digit == 4 else F(13915537, 10000000)
        require(record['digit'] == digit, 'Wrong pressure digit')
        require(record['N'] == expected_n, 'Wrong pressure grid')
        require(F(record['p']) == F(3, 2), 'Wrong moment order')
        require(F(record['pressure_bound']) == expected_g, 'Wrong pressure exponent')
        require(F(record['lambda']) == expected_lambda, 'Wrong upper eigenvalue')
        require(record['vector_denominator'] == pressure.V_DEN == 1 << 40,
                'Wrong vector denominator')
        require(record['weight_denominator'] == pressure.W_DEN == 1 << 40,
                'Wrong weight denominator')
        pressure.verify_record(record)

    print('PART 3: two-sided affine Hellinger profile certificates', flush=True)
    data = json.loads((HERE / 'hellinger_certificates.json').read_text())
    require(data['status'] == 'complete', 'Incomplete Hellinger certificate')
    require(data['mesh'] == hellinger.MESH == 200, 'Wrong Hellinger grid')
    require(data['fixed_denominator'] == hellinger.DEN == 1 << 28,
            'Wrong Hellinger denominator')
    require([row['digit'] for row in data['digits']] == list(range(5)),
            'Missing or repeated Hellinger digit')
    expected_targets = [F(19, 20), F(967, 1000), F(239, 250), F(481, 500), F(241, 250)]
    expected_deltas = [F(7, 1000), F(3, 625), F(3, 500), F(11, 2000), F(53, 10000)]
    require([hellinger.TARGETS[i] for i in range(5)] == expected_targets,
            'Hellinger target constants changed')
    require([hellinger.DELTAS[i] for i in range(5)] == expected_deltas,
            'Hellinger contraction exponents changed')
    for saved in data['digits']:
        current = hellinger.compute(saved['digit'], 200)
        require(current == saved, 'Recomputed Hellinger certificate differs')
        digit = saved['digit']
        require(F(current['certified_delta']) >=
                (F(53, 10000) if digit == 4 else F(1, 250)),
                'Contraction insufficient for the Type II input pair')
        print('VERIFIED Hellinger digits', digit, 'and', 9-digit,
              'K2 <=', current['K2_target'],
              'delta >=', current['certified_delta'], flush=True)

    print('PART 4: rebuild all 188 rational budget systems', flush=True)
    budgets.main()
    require(not any(name.split('.')[0] in {'numpy', 'scipy', 'sympy'}
                    for name in sys.modules), 'A numerical optimizer was imported')

    print('VERIFIED: 15 moment bounds, 5 Hellinger profiles, 188 rational systems.', flush=True)
    print('Elapsed seconds:', round(time.monotonic()-started, 3), flush=True)



if __name__ == '__main__':
    main()
