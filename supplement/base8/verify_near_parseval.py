"""Verify the shifted octal moment of order 127/64."""

from functools import lru_cache
from fractions import Fraction
import json
from math import isqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
from exact_trig import SCALE, ceil_div, cos_interval

NODAL_DATA = HERE / 'near_parseval_nodes.json'
TNUM, TDEN = 11481, 10000  # 1.1481
POWER, ROOT_LOG = 127, 7
COARSE, VERIFY, GSCALE = 256, 32768, 10000000


def verify():
    data = json.loads(NODAL_DATA.read_text())
    assert len(data) == 4 and all(len(v) == COARSE for v in data)
    R = VERIFY // COARSE
    assert R * COARSE == VERIFY

    @lru_cache(maxsize=None)
    def cos_hi(r):
        return cos_interval(r, 4 * VERIFY)[1]

    result = []
    for d, values in enumerate(data):
        assert all(isinstance(x, int) and x > 0 for x in values)
        nodal = [(R-b)*values[a]+b*values[(a+1)%COARSE]
                 for a in range(COARSE) for b in range(R)]
        G = GSCALE * R
        q = 8 * VERIFY
        # Coefficients in |sum_{a neq d} exp(2 pi i a t)|^2.
        diffs = [8-h-(d>=h)-(d+h<8) for h in range(1, 8)]
        assert all(a >= 0 for a in diffs)
        weights_hi = []
        for a in range(q):
            hi = 7 * SCALE
            for h, c in enumerate(diffs, 1):
                r = (h * a) % q
                hi += 2*c*cos_hi(min(r, q-r))
            u = min(SCALE, ceil_div(hi, 49))
            assert u >= 0
            power = SCALE * u ** POWER
            y = power
            for _ in range(ROOT_LOG):
                y = isqrt(y)
            weights_hi.append(y + (y ** (1 << ROOT_LOG) < power))

        max_delta = None
        for i, gi in enumerate(nodal):
            up = 0
            for j in range(8):
                ell = i + j * VERIFY
                a, b = divmod(ell, 8)
                up += weights_hi[ell] * ((8 - b) * nodal[a]
                                        + b * nodal[(a + 1) % VERIFY])
            # Common denominator: 8 * SCALE * G * TDEN.
            delta = TDEN * up - (8 * SCALE * TNUM) * gi
            max_delta = delta if max_delta is None else max(delta, max_delta)
        lip_g = Fraction(VERIFY * max(abs(nodal[a] - nodal[(a+1)%VERIFY])
                                     for a in range(VERIFY)), G)
        # |phi| <= 1. Center the phases around 7/2; its Lipschitz
        # constant is <= (2*pi/7)*16 < 15, so |phi|^(127/64) is
        # Lipschitz with constant <=30. Apply the bound to g-1.
        # The exact square operator has T_2 1=8/7. The derivative
        # of x^(127/64)-x^2 is bounded by 2/64 for 0<=x<=1;
        # hence T_s 1 has Lipschitz constant <=30/64.
        lam = Fraction(TNUM, TDEN)
        sup_dev = Fraction(max(abs(v-G) for v in nodal), G)
        guard = (Fraction(30,64) + 30*sup_dev + (1+lam)*lip_g) / (2*VERIFY)
        margin = -Fraction(max_delta, 8 * SCALE * G * TDEN) - guard
        assert margin > 0, (d, str(margin))
        result.append({'digit': d, 'continuous_margin': str(margin)})

    # lambda < (8/7)*8^(0.141 / 64).
    assert (7*TNUM)**64000 < 8**141 * (8*TDEN)**64000
    print(json.dumps({'lambda': '11481/10000', 'verified': result,
                      'exponent_increment': '141/64000'}, indent=2))



if __name__ == '__main__':
    verify()
