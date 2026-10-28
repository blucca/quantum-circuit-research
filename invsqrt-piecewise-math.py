#!/usr/bin/env python3
"""Integer-Horner accuracy probe for a direct piecewise-cubic inverse root.

Scope: the public 22-bit arithmetic contract, before reversible circuit synthesis.
The fixed choices are 32 mantissa intervals, degree 3, 28 coefficient bits.
Horner multiplies use floor shifts, including for negative intermediate values.
No competition dataset, circuit or secret-salted verifier is used here.
"""
import json
import os
from pathlib import Path
import resource
import time

os.environ['OPENBLAS_NUM_THREADS'] = '1'
resource.setrlimit(resource.RLIMIT_CPU, (30, 30))
resource.setrlimit(resource.RLIMIT_AS, (1024 * 1024**2, 1024 * 1024**2))
import numpy as np
from numpy.polynomial import Chebyshev, Polynomial

MAX_X, SEGMENTS, FRAC, MANTISSA_BITS = 3_133_452, 32, 28, 21
CHUNK = 1 << 18
OUT = Path(os.environ.get('PROBE_OUT', 'temp/invsqrt-probe'))
OUT.mkdir(parents=True, exist_ok=True)


def make_coefficients():
    coefficients = np.empty((2, SEGMENTS, 4), dtype=np.int64)
    h = 1.0 / SEGMENTS
    for i in range(SEGMENTS):
        left = 1 + i * h
        fit = Chebyshev.interpolate(lambda z: 1 / np.sqrt(left + (z + 1) * h / 2), 3)
        local = fit.convert(kind=Polynomial)(Polynomial([-1, 2 / h])).coef
        for odd in (0, 1):
            coefficients[odd, i] = np.rint(local * (2.0**FRAC) / np.sqrt(2.0**odd)).astype(np.int64)
    assert coefficients.shape == (2, 32, 4)
    assert np.abs(coefficients).max() < 2**31
    return coefficients


def scalar_horner(x, coefficients):
    k = int(x).bit_length() - 1
    m = int(x) << (MANTISSA_BITS - k)
    segment = (m - (1 << MANTISSA_BITS)) >> (MANTISSA_BITS - 5)
    delta = m & ((1 << (MANTISSA_BITS - 5)) - 1)
    c = [int(v) for v in coefficients[k & 1, segment]]
    acc = c[3]
    for j in (2, 1, 0):
        acc = ((acc * delta) >> MANTISSA_BITS) + c[j]
    return acc, k // 2


def main():
    cpu0 = time.process_time()
    coefficients = make_coefficients()
    (OUT / 'coefficients.json').write_text(json.dumps(coefficients.tolist()))
    maximum_relative = 0.0
    maximum_ratio = 0.0
    worst_x = None
    minimum_phase_slack = float('inf')
    max_abs_product, max_abs_accumulator = 0, 0
    checked = 0
    for start in range(1, MAX_X + 1, CHUNK):
        x = np.arange(start, min(MAX_X + 1, start + CHUNK), dtype=np.int64)
        # log2 is exact at powers of two; all inputs are <2^22 and stored exactly.
        k = np.floor(np.log2(x)).astype(np.int64)
        assert np.all((1 << k) <= x) and np.all(x < (1 << (k + 1)))
        mantissa = np.left_shift(x, MANTISSA_BITS - k)
        index = np.right_shift(mantissa - (1 << MANTISSA_BITS), MANTISSA_BITS - 5)
        delta = mantissa & ((1 << (MANTISSA_BITS - 5)) - 1)
        c = coefficients[k & 1, index]
        acc = c[:, 3].copy()
        for j in (2, 1, 0):
            product = acc * delta
            max_abs_product = max(max_abs_product, int(np.abs(product).max()))
            acc = np.right_shift(product, MANTISSA_BITS) + c[:, j]
            max_abs_accumulator = max(max_abs_accumulator, int(np.abs(acc).max()))
        for j in (0, len(x) // 2, len(x) - 1):
            value, shift = scalar_horner(int(x[j]), coefficients)
            assert value == int(acc[j]) and shift == int(k[j] // 2)
        approximation = np.ldexp(acc.astype(np.float64), -FRAC - k // 2)
        exact = 1 / np.sqrt(x.astype(np.float64))
        absolute = np.abs(approximation - exact)
        relative = absolute / exact
        pos = int(relative.argmax())
        if float(relative[pos]) > maximum_relative:
            maximum_relative, worst_x = float(relative[pos]), int(x[pos])
        maximum_ratio = max(maximum_ratio, float(relative[pos]) * 2**24)
        minimum_phase_slack = min(minimum_phase_slack, float((exact * 2**-24 - absolute).min()))
        checked += len(x)
        state = {'phase': 'running', 'checked': checked, 'max_relative_error': maximum_relative,
                 'max_tolerance_ratio': maximum_ratio, 'worst_x': worst_x,
                 'minimum_phase_slack_turns': minimum_phase_slack}
        (OUT / 'state.json').write_text(json.dumps(state, indent=2))
    assert checked == MAX_X and max_abs_product < 2**63
    # Check the published arithmetic tolerance; report headroom separately.
    assert maximum_ratio < 1.0 and minimum_phase_slack > 0
    state.update(phase='complete', passed=True, eightfold_headroom_passed=maximum_ratio < 0.125,
                 tolerance_headroom=1 / maximum_ratio, max_abs_product=max_abs_product,
                 max_abs_accumulator=max_abs_accumulator, cpu_seconds=time.process_time() - cpu0,
                 peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                 algorithm={'segments': SEGMENTS, 'degree': 3, 'coefficient_fraction_bits': FRAC,
                            'mantissa_fraction_bits': MANTISSA_BITS, 'horner_rounding': 'floor'},
                 circuit_status='arithmetic only; reversible synthesis and official verification are next')
    (OUT / 'state.json').write_text(json.dumps(state, indent=2))
    print(json.dumps(state, indent=2))


if __name__ == '__main__':
    main()
