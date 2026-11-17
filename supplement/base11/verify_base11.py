#!/usr/bin/env python3
"""Verify every base-eleven finite input, using only the standard library."""
from fractions import Fraction as F
from pathlib import Path
import runpy

assert __debug__, 'Run without -O.'
HERE=Path(__file__).resolve().parent

def check_parameters():
    assert 10**500 > 11**477           # c_11 > .954
    assert 10**100 < 11**99            # c_11 < .99
    assert 10**21 > 11**20             # c_11/2 > 10/21
    assert 1001**5000 < 11**3*1000**5000
    gamma=F(50,77); c=F(477,500)
    assert gamma*c/(c+gamma*(1-c)) > F(629,1000)
    assert F(27,77)+F(10,21) < F(5,6)
    assert 1+F(10,21) < F(3,2)
    assert F(23,80)*(2-F(235,154))-F(59,433) == F(127,5334560)
    assert 24**250*11**3 <= 25**250
    print('VERIFIED all scalar base-11 inequalities',flush=True)

if __name__=='__main__':
    check_parameters()
    runpy.run_path(str(HERE/'verify_moments.py'),run_name='__main__')
    runpy.run_path(str(HERE/'verify_hellinger.py'),run_name='__main__')
