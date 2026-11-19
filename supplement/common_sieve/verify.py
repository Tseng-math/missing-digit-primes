"""Compile one directed simplex enumerator and verify the requested windows."""
from fractions import Fraction as Q
from math import factorial
from os import environ
from pathlib import Path
from subprocess import run
import sys
from tempfile import TemporaryDirectory


SCALE = 10**15
WINDOWS = {
    'six': (4230, 4800, (
        (1000, 3, (776428472603853, 170088578625634)),
        (400, 5, (36207776889770, 2604877435452)),
        (200, 7, (236130635415, 7040715097, 133961892)),
        (160, 10, (11398446, 113413, 55354, 59498, 110569, 133507, 139568, 45760)),
    )),
    'seven_pairs': (4200, 4800, (
        (640, 3, (740985917419842, 155577002095249, 24073031839013, 1371042713357,
                  56158753041, 1189984604, 31804735)),
    )),
    'nine': (3155, 3691, (
        (512, 3, (676842122671633, 98795867966415, 53380320592348, 5037807356947,
                  131073418223, 2250013292, 51176506, 9579102, 4896745, 5262174)),
    )),
}


def verify(binary, name):
    left, right, batches = WINDOWS[name]
    eta = Q(right-left, 10000)
    ranks = []
    total = 0
    for mesh, first, expected in batches:
        last = first + len(expected) - 1
        args = [left, right, mesh, first, last]
        output = run([binary, *map(str, args)], check=True, capture_output=True, text=True).stdout
        lines = output.splitlines()
        assert len(lines) == len(expected) + 1
        for rank, bound, line in zip(range(first, last+1), expected, lines):
            assert line.split()[:2] == [f'rank={rank}', f'upper={bound}'], (name, line, bound)
            ranks.append(rank)
        assert lines[-1] == f'total={sum(expected)}'
        total += sum(expected)
        print(f'{name} mesh={mesh}\n{output}', end='', flush=True)
    assert ranks == list(range(3, last+1))
    tail = Q(0)
    for rank in range(last+1, 1 + 1//eta):
        volume = (1-rank*eta)**(rank-1) / (factorial(rank)*factorial(rank-1))
        product = eta**rank if name == 'nine' else eta**(rank-1)*(1-(rank-1)*eta)
        tail += 2**(rank-2) * volume / product
    loss = Q(total, SCALE) + tail
    if name == 'six':
        assert last == 1//eta == 17 and not tail and total == 985573022923228
        assert Q(3,25)*(1-loss) == Q('.00173123724921264') > Q(1,600)
        cap = '0.985573022923228 (inclusive)'
    elif name == 'seven_pairs':
        assert last == 9 and 1//eta == 16 and total == 922064374609841
        assert tail < Q(11298778668719, SCALE)
        limit = Q(933363153278560, SCALE)
        assert limit - Q(1,SCALE) <= loss < limit < Q('.934')
        assert Q(35,216)*(1-limit) > Q('.01')
        cap = '0.933363153278560 (strict)'
    else:
        assert last == 12 and 1//eta == 18 and total == 834189512933385
        assert tail < Q('.000013819') and loss < Q('.834204')
        assert Q('.165796')*Q(15,16) > Q(2,13) and Q(2,13)/8 == Q(1,52)
        assert Q(9,8)/eta < 21 and 22*Q(8,7)*8 < 400
        cap = '0.834204 (strict)'
    print(f'PASS {name}: exact tail={tail}; loss bound {cap}', flush=True)


if __name__ == '__main__':
    names = sys.argv[1:] or list(WINDOWS)
    if not __debug__ or not set(names) <= WINDOWS.keys() or len(set(names)) != len(names):
        raise SystemExit('Enable assertions and select distinct known windows.')
    root = Path(__file__).resolve().parent
    with TemporaryDirectory(dir=root) as directory:
        binary = Path(directory)/'simplex'
        env = dict(environ, TMPDIR=directory, TMP=directory, TEMP=directory, CLANG_MODULE_CACHE_PATH=directory)
        run([env.get('CXX','c++'), '-std=c++17', '-O3', str(root/'simplex.cpp'), '-o', binary],
            cwd=directory, env=env, check=True)
        for name in names:
            verify(binary, name)
