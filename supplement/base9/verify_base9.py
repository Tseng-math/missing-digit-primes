#!/usr/bin/env python3
"""Recompute the finite certificates used for all nine forbidden digits.

Only Python's standard library and a C++17 compiler supporting unsigned
__int128 are needed.  Positive Fourier supersolutions, scalar inequalities,
and every N=512 complete-tuple cell are checked afresh.  Generated reports
and the compiled executable are confined to an automatically removed
temporary directory; no stored verification output is used as evidence.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
SOURCES = (
    'verify_moments.py',
    'moment_certificates.json',
    'certify_full_tuple_sieve.cpp',
)
EXPECTED_TUPLE_NUMERATORS = {
    3: 676842122666188091,
    4: 98795867909592915,
    5: 129047629062754398,
    6: 18176369591635255,
    7: 535666342176341,
    8: 7612230320382,
    9: 76343360570,
    10: 529406458,
    11: 7238880,
    12: 3166020,
    13: 1492040,
    14: 453596,
    15: 252616,
    16: 115570,
    17: 0,
    18: 0,
}


def check_parameters():
    sigma, s, g_s = F(127, 200), F(3, 2), F(6307, 40000)
    kappa, g4 = g_s / (2-s), ((1-sigma)+2*g_s)/3
    c0, y = F(473, 500), F(633, 1000)
    delta0, omega = F(1, 2000), F(1, 10**6)
    a0, b0, u0 = F(6309, 20000), F(7383, 20000), F(1477, 4000)
    eps0, xi = F(1, 10**6), F(1, 10**5)
    a, b, eta = F(631, 2000), F(3691, 10000), F(67, 1250)
    assert 8**500 > 9**473
    assert kappa == F(6307, 20000) and g4 == F(13607, 60000)
    assert F(4, 3) <= s < 2 and kappa < F(1, 3)
    assert 2*kappa < y < min(F(2, 3), sigma, 1-F(3, 2)*g4)
    rank_margins = {
        'first_energy': sigma+c0/2-1,
        'middle_energy': 1-y-3*g4/2,
        'third_energy': 2*sigma-2*y,
        'determinant': 2-3*y,
    }
    assert rank_margins == {
        'first_energy': F(27, 250),
        'middle_energy': F(1073, 40000),
        'third_energy': F(1, 250),
        'determinant': F(101, 1000),
    }
    assert delta0 < min(rank_margins.values())/4
    rho, tau = 2*sigma-1, c0+2*sigma-2
    assert rho-omega > F(1, 5) and tau-omega > F(1, 5)
    assert omega < min(rho/2, tau/2, delta0*(1-y)/(2*y))
    global_saving = delta0*(1-y)-y*omega
    assert global_saving == F(182867, 10**9) > F(1, 6000)
    typeii_margins = {
        'generic': 1-u0-2*kappa-2*(s-1)*eps0,
        'K_prime_grows': u0-b0,
        'rank_two_range': y-(1-u0),
        'nonzero_line_count': 2*a0+u0-1,
        'line_first': 2*a0-2*kappa-(3*s-4)*eps0-xi,
        'line_second': 3*a0+1-6*kappa-6*(s-1)*eps0-2*xi,
        'tail': (2-s)*eps0/2,
    }
    assert typeii_margins == {
        'generic': F(49, 10**6),
        'K_prime_grows': F(1, 10000),
        'rank_two_range': F(9, 4000),
        'nonzero_line_count': F(3, 20000),
        'line_first': F(379, 2000000),
        'line_second': F(54227, 10**6),
        'tail': F(1, 4000000),
    }
    assert min(typeii_margins.values()) > 0
    assert kappa < a0 < a < b < b0 < 1-2*kappa
    assert eta == b-a and int(1/eta) == 18
    assert F(1, 4) < a < F(1, 3) < b < F(1, 2)
    assert 2*a+b > 1 and F(5, 2)*a < 1
    assert 1-b+F(1, 500) < sigma < c0
    print('PASS: all exact base-9 Fourier, geometry, Type II, and sieve parameters', flush=True)
    return {
        'rank_margins': {key: str(value) for key, value in rank_margins.items()},
        'typeii_margins': {key: str(value) for key, value in typeii_margins.items()},
        'global_rank_saving': str(global_saving),
    }


def check_moment_report(path, certificate):
    report = json.loads(path.read_text())
    assert report['status'] == 'PASS: exact integer/rational continuous supersolution certificates'
    assert report['base'] == 9 and report['forbidden_digits'] == list(range(9))
    assert report['moment_one_discrete_exponent'] == '73/200'
    assert report['moment_one_continuum_decay'] == '127/200'
    assert report['moment_three_halves_discrete_exponent'] == '6307/40000'
    assert report['moment_three_halves_continuum_decay'] == '33693/40000'
    assert report['kappa'] == '6307/20000'
    assert report['interpolated_four_thirds_discrete_exponent'] == '13607/60000'
    assert report['certificate_sha256'] == hashlib.sha256(certificate.read_bytes()).hexdigest()
    rows = report['results']
    assert len(rows) == 10
    assert {(row['forbidden_digit'], F(row['moment'])) for row in rows} == {
        (digit, moment) for digit in range(5) for moment in (F(1), F(3, 2))
    }
    for row in rows:
        assert row['reflected_digit'] == 8-row['forbidden_digit']
        assert F(row['minimum_g']) > 0 and F(row['strict_margin']) > 0
        assert F(row['lambda_bound']) == (
            F(2477, 10000) if F(row['moment']) == 1 else F(314229, 2000000)
        )
    print('PASS: all 10 continuous supersolutions, covering all 9 digits by reflection', flush=True)


def check_tuple_report(path):
    report = json.loads(path.read_text())
    assert report['N'] == 512 and int(report['denominator']) == 10**18
    assert (report['a'], report['b'], report['eta']) == ('631/2000', '3691/10000', '67/1250')
    rows = report['results']
    assert [row['k'] for row in rows] == list(range(3, 19))
    current = {row['k']: int(row['upper_numerator']) for row in rows}
    assert current == EXPECTED_TUPLE_NUMERATORS
    numerator = sum(current.values())
    assert numerator == int(report['total_upper_numerator']) == 923405344688153132
    total = F(numerator, 10**18)
    assert total == F(report['total_upper_decimal']) < F(231, 250)
    assert 1-total > F(19, 250)
    assert F(15, 16)*(1-total) > F(1, 40)
    assert F(1, 40)/8 == F(1, 320)
    assert F(9, 8)/F(67, 1250) < 21 < 22
    assert F(22)*F(8, 7)*8 < 400
    print('PASS: all 16 complete-tuple dimensions r=3,...,18 at N=512', flush=True)
    print('PASS: total numerator 923405344688153132 / 1000000000000000000 < 231/250', flush=True)
    return numerator


def main():
    if sys.flags.optimize:
        raise RuntimeError('Run without -O: all exact certificate assertions must remain enabled.')
    compiler = shutil.which('g++')
    if compiler is None:
        raise RuntimeError('g++ with C++17 and unsigned __int128 support is required.')
    parameters = check_parameters()
    with tempfile.TemporaryDirectory(prefix='verification-', dir=HERE) as td:
        work = Path(td)
        for name in SOURCES:
            shutil.copyfile(HERE/name, work/name)
        subprocess.run([sys.executable, '-B', '-S', str(work/SOURCES[0])],
                       cwd=work, check=True)
        check_moment_report(work/'moment_verification.json', work/SOURCES[1])
        executable = work/'certify_full_tuple_sieve'
        subprocess.run([compiler, '-O3', '-std=c++17', str(work/SOURCES[2]),
                        '-o', str(executable)], cwd=work, check=True,
                       env={**os.environ, 'TMPDIR': td, 'TMP': td, 'TEMP': td})
        output = work/'full_tuple_sieve_n512.json'
        subprocess.run([str(executable), '--n', '512', '--kmin', '3', '--kmax', '18',
                        '--out', str(output)], cwd=work, check=True)
        numerator = check_tuple_report(output)
    print(json.dumps({
        'status': 'PASS',
        'finite_premises': '10 continuous Fourier supersolutions and all 16 complete-tuple upper sums',
        'tuple_loss_numerator': numerator,
        'tuple_loss_denominator': 10**18,
        'strict_loss_bound': '231/250',
        **parameters,
    }, indent=2))


if __name__ == '__main__':
    main()
