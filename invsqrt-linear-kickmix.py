#!/usr/bin/env python3
"""Frozen 1024-cell biased-linear circuit with one exact signed product.

Joint QROM holds unsigned29 intercept and signed19 local cell change.
Product is signed30, extended to 40 before adding the intercept.  All bits
survive until their exact inverse; the existing public HMR adders clean carries.
"""
from fractions import Fraction
import importlib.util
import json
import math
import os
from pathlib import Path
import time

import kickmix as km
import numpy as np

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('cubic_helpers', ROOT/'invsqrt-kickmix.py')
helpers = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helpers)
arr, highest_bit = helpers.arr, helpers.highest_bit
COEFF = json.loads((ROOT/'results/invsqrt-linear-coefficients.json').read_text())
MATH = json.loads((ROOT/'results/invsqrt-linear-math.json').read_text())
MAX_X = 3_133_452
OUT = Path(os.environ.get('PROBE_OUT', 'temp/invsqrt-linear-circuit-probe'))


def signed_product(b, d, delta, p, *, inverse=False):
    assert (len(d), len(delta), len(p)) == (19, 11, 40)
    if inverse:
        b.cx(p[29], p[30:])
        b.iadd(delta, target=p[19:30], control=d[18])
        for i in range(10, -1, -1):
            b.isub(d, target=p[i:30], control=delta[i])
    else:
        for i in range(11):
            b.iadd(d, target=p[i:30], control=delta[i])
        b.isub(delta, target=p[19:30], control=d[18])
        b.cx(p[29], p[30:])


def make_circuit():
    b = km.CircuitBuilder()
    b.set_max_qubits(150)
    x = b.create_quantum_register(22, name='x')
    m, k = b.alloc_qubits(22), b.alloc_qubits(5)
    highest_bit(b, x, k)
    b.cx(x, m)
    b.left_rotate(m, 21)
    b.right_rotate(m, k)
    delta = m[:11]
    address = arr(list(m[11:21])+[k[0]])
    table = [int(COEFF[o][s][0]) | ((int(COEFF[o][s][1]) % (1 << 19)) << 29)
             for o in range(2) for s in range(1024)]
    coefficients = b.alloc_qubits(48)
    b.init_lookup(table=table, address=address, target=coefficients)
    a, d = coefficients[:29], coefficients[29:]
    p = b.alloc_qubits(40)
    signed_product(b, d, delta, p)
    b.iadd(a, target=p[11:])
    zeros = b.alloc_qubits(10)
    scaled = arr(list(zeros)+list(p))
    b.right_rotate(scaled, k[1:])
    # z_pow is in half-turns; the most significant wire contributes half a turn.
    assert len(scaled) == 50
    for i in range(50):
        b.z_pow(scaled[i], Fraction(1, 1 << (49-i)))
    b.left_rotate(scaled, k[1:])
    b.free(zeros)
    b.isub(a, target=p[11:])
    signed_product(b, d, delta, p, inverse=True)
    b.free(p)
    b.del_lookup(table=table, address=address, target=coefficients)
    b.free(coefficients)
    b.left_rotate(m, k)
    b.right_rotate(m, 21)
    b.cx(x, m)
    b.free(m)
    highest_bit(b, x, k)
    b.free(k)
    return b.finish_circuit()


def scalar(x):
    k = int(x).bit_length()-1
    m = int(x) << (21-k)
    s, delta = (m >> 11) & 1023, m & 2047
    a, d = map(int, COEFF[k & 1][s])
    return Fraction((a << 11)+d*delta, 1 << (40+k//2))


def probe(circuit):
    xs = {1, 2, 3, MAX_X, MATH['worst_x']}
    for k in range(22):
        for s in (0, 1, 511, 512, 1023, 1024):
            boundary = ((1024+s) << k)//1024
            for x in (boundary-1, boundary, boundary+1):
                if 1 <= x <= MAX_X:
                    xs.add(x)
    xs.update(int(x) for x in np.random.default_rng(262).integers(1, MAX_X+1, 256))
    xs = sorted(xs)
    shots = 64
    sim = km.Simulator(batch_size=shots, counting_operations=True)
    sim.use_same_registers_as(circuit)
    scratch = arr(km.q(q) for q in range(22, circuit.num_qubits))
    for checked, x in enumerate(xs, 1):
        sim.clear_for_shot()
        bits = ((x >> np.arange(22)) & 1).astype(bool)[:, None].repeat(shots, axis=1)
        sim.write_across_shots('x', bits)
        sim.do(circuit)
        assert np.array_equal(sim.read_across_shots('x'), bits), ('input changed', x)
        assert not sim.read_across_shots(scratch).any(), ('dirty scratch', x)
        expected = scalar(x) % 1
        for shot in range(shots):
            assert Fraction(sim.read_shot_phase(shot), 2) % 1 == expected, ('phase mismatch', x, shot)
        if checked % 64 == 0:
            (OUT/'state.json').write_text(json.dumps({'phase': 'testing', 'inputs': checked, 'shots': checked*shots}))
    stats = circuit.operation_statistics()
    sampled = sim.operation_statistics()
    fine = stats.num_z_pow(min_clifford_hierarchy_level=4)
    slack = MATH['minimum_phase_slack_turns']
    eps = math.pi*slack/fine
    rotation_t = fine*(0.56*math.log2(1/eps)+5.3)
    n = len(xs)*shots
    mean_toffoli, mean_t = sampled.num_toffoli_like/n, sampled.num_t_like/n
    cost = mean_toffoli+(mean_t+rotation_t)/4
    volume = cost*circuit.num_qubits
    return {'phase': 'complete', 'sampled_circuit_passed': True,
            'inputs': len(xs), 'shots': n, 'input_restored': True, 'scratch_clean': True,
            'phase_check': 'exact integer linear numerator on every sampled HMR shot',
            'qubits': circuit.num_qubits, 'sample_mean_toffoli_like': mean_toffoli,
            'sample_mean_t_like': mean_t, 'non_clifford_rotations': fine,
            'rotation_epsilon': eps, 'rotation_t_gates_paper_formula': rotation_t,
            'toffoli_equivalents_paper_formula': cost, 'qubits_times_cost': volume,
            'reference_qubits_times_cost': 260*7988, 'fraction_of_reference': volume/(260*7988),
            'resource_gate_passed': volume <= 0.5*260*7988,
            'official_verification': 'not run; official generator and full 200540928 shots are a later gate',
            'operation_counts': {key: int(value) for key, value in stats.items()}}


if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.process_time()
    (OUT/'state.json').write_text(json.dumps({'phase': 'generating'}))
    try:
        circuit = make_circuit()
        circuit.to_file(OUT/'invsqrt.kmx')
        result = probe(circuit)
        result['cpu_seconds'] = time.process_time()-t0
        (OUT/'state.json').write_text(json.dumps(result, indent=2))
        print(json.dumps(result, indent=2), flush=True)
    except Exception as e:
        (OUT/'failure.json').write_text(json.dumps({'phase': 'failed', 'error': repr(e)}))
        raise
