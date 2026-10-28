#!/usr/bin/env python3
"""One frozen biased-linear arithmetic witness; no table/width search.

1024 cells, 29 fractional bits in both local coefficients.  The local cell
change fits signed19, and the exact numerator fits unsigned40.  The bias is
half the curvature upper bound on the chord error in each cell.
"""
import json
import os
from pathlib import Path
import resource
import time

os.environ['OPENBLAS_NUM_THREADS'] = '1'
resource.setrlimit(resource.RLIMIT_CPU, (30, 30))
resource.setrlimit(resource.RLIMIT_AS, (768 * 1024**2, 768 * 1024**2))
import numpy as np

MAX_X, CELLS, FRAC, LOW_BITS = 3_133_452, 1024, 29, 11
CHUNK = 1 << 18
OUT = Path(os.environ.get('PROBE_OUT', 'temp/invsqrt-linear-math-probe'))
OUT.mkdir(parents=True, exist_ok=True)


def coefficients():
    a = 1 + np.arange(CELLS, dtype=np.float64) / CELLS
    h = 1 / CELLS
    f = 1 / np.sqrt(a)
    bias = (3 / 4) * a**(-2.5) * h*h / 16
    change = 1 / np.sqrt(a+h) - f
    table = np.rint(np.stack([
        np.stack([f-bias, change], axis=-1) / np.sqrt(2**odd)
        for odd in (0, 1)]) * 2**FRAC).astype(np.int64)
    assert np.all((0 <= table[..., 0]) & (table[..., 0] < 2**29))
    assert np.all((-2**18 <= table[..., 1]) & (table[..., 1] < 2**18))
    return table


def scalar(x, table):
    k = int(x).bit_length()-1
    m = int(x) << (21-k)
    s, delta = (m >> LOW_BITS) & (CELLS-1), m & ((1 << LOW_BITS)-1)
    a, d = map(int, table[k & 1, s])
    numerator = (a << LOW_BITS) + d*delta
    assert 0 <= numerator < 2**40
    return numerator, 40+k//2


def main():
    t0 = time.process_time()
    table = coefficients()
    (OUT/'coefficients.json').write_text(json.dumps(table.tolist()))
    state = {'phase': 'running', 'algorithm': {'cells': CELLS, 'frac': FRAC,
             'intercept_bits': 29, 'local_change_signed_bits': 19,
             'delta_bits': 11, 'numerator_bits': 40, 'rounding': 'nearest coefficients; exact product'}}
    maximum = 0.0
    minimum_slack = float('inf')
    min_n, max_n, max_product = 2**40, 0, 0
    for start in range(1, MAX_X+1, CHUNK):
        x = np.arange(start, min(MAX_X+1, start+CHUNK), dtype=np.int64)
        k = np.floor(np.log2(x)).astype(np.int64)
        assert np.all((1 << k) <= x) and np.all(x < (1 << (k+1)))
        m = np.left_shift(x, 21-k)
        s, delta = (m >> LOW_BITS) & (CELLS-1), m & ((1 << LOW_BITS)-1)
        c = table[k & 1, s]
        product = c[:, 1]*delta
        numerator = (c[:, 0] << LOW_BITS)+product
        assert np.all((0 <= numerator) & (numerator < 2**40))
        assert np.all((-2**29 <= product) & (product < 2**29))
        for i in (0, len(x)//2, len(x)-1):
            assert scalar(int(x[i]), table) == (int(numerator[i]), int(40+k[i]//2))
        approximation = np.ldexp(numerator.astype(np.float64), -40-k//2)
        exact = 1/np.sqrt(x.astype(np.float64))
        absolute = np.abs(approximation-exact)
        relative = absolute/exact
        i = int(relative.argmax())
        if float(relative[i]) > maximum:
            maximum, worst_x = float(relative[i]), int(x[i])
        minimum_slack = min(minimum_slack, float((exact*2**-24-absolute).min()))
        min_n, max_n = min(min_n, int(numerator.min())), max(max_n, int(numerator.max()))
        max_product = max(max_product, int(np.abs(product).max()))
        state.update(checked=int(x[-1]), maximum_relative_error=maximum, worst_x=worst_x,
                     minimum_phase_slack_turns=minimum_slack)
        (OUT/'state.json').write_text(json.dumps(state, indent=2))
    passed = maximum < 2**-24 and minimum_slack > 0
    state.update(phase='complete', passed=passed, max_tolerance_ratio=maximum*2**24,
                 numerator_range=[min_n, max_n], max_absolute_product=max_product,
                 cpu_seconds=time.process_time()-t0,
                 peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                 circuit_status='arithmetic witness only')
    (OUT/'state.json').write_text(json.dumps(state, indent=2))
    print(json.dumps(state, indent=2))
    assert passed, 'Frozen linear integer arithmetic failed the actual contract'


if __name__ == '__main__':
    main()
