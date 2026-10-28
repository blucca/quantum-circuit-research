#!/usr/bin/env python3
"""Exact reversible implementation of the already measured 28-bit Horner probe.

Uses public kickmix 0.1.dev0. Keeps the three product histories for transparent
cleanup. This is a first measured circuit, not an optimized competition entry.
All multiplication is exact; no carry truncation or probabilistic arithmetic.
"""
from fractions import Fraction
import json
import math
import os
from pathlib import Path
import time

import kickmix as km
import numpy as np

ROOT = Path(__file__).resolve().parent
COEFF = json.loads((ROOT / 'results/invsqrt-piecewise-coefficients.json').read_text())
MAX_X = 3_133_452
OUT = Path(os.environ.get('PROBE_OUT', 'temp/invsqrt-circuit-probe'))


def arr(xs):
    return km.array(list(xs))


def highest_bit(b, x, k):
    """XOR floor(log2(x)) into k, x!=0. Prefix work is exactly cleared."""
    prefix = []
    above = False
    for j in range(21, 0, -1):
        out = b.alloc_qubits(1)[0]
        b.cx(x[j], out)
        if above is not False:
            b.cx(above, out)
            b.ccx(x[j], above, out)
        prefix.append((x[j], above, out))
        above = out
    # prefix[20-j] is OR(x[j+1:]) for j<21.
    for j in range(1, 22):
        prev = False if j == 21 else prefix[20-j][2]
        if prev is False:
            for bit in range(5):
                if j >> bit & 1:
                    b.cx(x[j], k[bit])
        else:
            match = b.alloc_qubits(1)[0]
            b.x(prev)
            b.ccx(x[j], prev, match)
            for bit in range(5):
                if j >> bit & 1:
                    b.cx(match, k[bit])
            b.del_and(x[j], prev, match)
            b.x(prev)
            b.free(arr([match]))
    for a, prev, out in reversed(prefix):
        if prev is not False:
            b.ccx(a, prev, out)
            b.cx(prev, out)
        b.cx(a, out)
        b.free(arr([out]))


def signed_product(b, a, delta, p, *, inverse=False):
    """p += signed29(a)*unsigned16(delta) from zero; inverse clears p.

    Compute the unsigned 45-bit product, apply the sign correction, then
    sign-extend by five distinct wires. Operands stay unchanged, and iadd's
    temporary carries use kickmix's measurement-based exact cleanup.
    """
    assert len(a) == 29 and len(delta) == 16 and len(p) == 50
    if inverse:
        b.cx(p[44], p[45:])
        b.iadd(delta, target=p[29:45], control=a[28])
        for i in range(15, -1, -1):
            b.isub(a, target=p[i:45], control=delta[i])
    else:
        for i in range(16):
            b.iadd(a, target=p[i:45], control=delta[i])
        b.isub(delta, target=p[29:45], control=a[28])
        b.cx(p[44], p[45:])


def make_circuit():
    b = km.CircuitBuilder()
    # A bounded first circuit; enough workspace for the library's clean adders.
    b.set_max_qubits(300)
    x = b.create_quantum_register(22, name='x')
    m = b.alloc_qubits(22)
    k = b.alloc_qubits(5)
    highest_bit(b, x, k)
    b.cx(x, m)
    b.left_rotate(m, 21)
    b.right_rotate(m, k)
    delta = m[:16]
    address = arr(list(m[16:21]) + [k[0]])
    tables = [[int(COEFF[o][s][j]) % (1 << 29)
               for o in range(2) for s in range(32)] for j in range(4)]
    seed = b.alloc_qubits(29)
    b.init_lookup(table=tables[3], address=address, target=seed)
    a = seed
    history = []
    for j in (2, 1, 0):
        p = b.alloc_qubits(50)
        signed_product(b, a, delta, p)
        coefficient = b.alloc_qubits(29)
        b.init_lookup(table=tables[j], address=address, target=coefficient)
        b.iadd(coefficient, target=p[21:])
        b.del_lookup(table=tables[j], address=address, target=coefficient)
        b.free(coefficient)
        history.append((j, a, p))
        a = p[21:]
    # A0 / 2^(28 + k//2). Ten clean low wires preserve every shifted bit.
    zeros = b.alloc_qubits(10)
    scaled = arr(list(zeros) + list(a))
    b.right_rotate(scaled, k[1:])
    for i in range(38):
        b.z_pow(scaled[i], Fraction(1, 1 << (37-i)))
    b.left_rotate(scaled, k[1:])
    b.free(zeros)
    for j, a, p in reversed(history):
        coefficient = b.alloc_qubits(29)
        b.init_lookup(table=tables[j], address=address, target=coefficient)
        b.isub(coefficient, target=p[21:])
        b.del_lookup(table=tables[j], address=address, target=coefficient)
        b.free(coefficient)
        signed_product(b, a, delta, p, inverse=True)
        b.free(p)
    b.del_lookup(table=tables[3], address=address, target=seed)
    b.free(seed)
    b.left_rotate(m, k)
    b.right_rotate(m, 21)
    b.cx(x, m)
    b.free(m)
    highest_bit(b, x, k)
    b.free(k)
    return b.finish_circuit()


def scalar_horner(x):
    k = x.bit_length()-1
    m = x << (21-k)
    s = (m-(1 << 21)) >> 16
    delta = m & 65535
    c = COEFF[k & 1][s]
    a = c[3]
    for j in (2,1,0):
        a = ((a*delta) >> 21) + c[j]
    return Fraction(a, 1 << (28+k//2))


def probe(circuit):
    # Edge/segment boundaries and random cases, 64 independent HMR shots each.
    xs = {1,2,3,MAX_X,999371}
    for k in range(22):
        for s in (0,1,15,16,30,31,32):
            boundary = ((32+s) << k) // 32
            for x in (boundary-1,boundary,boundary+1):
                if 1 <= x <= MAX_X:
                    xs.add(x)
    xs.update(int(x) for x in np.random.default_rng(260).integers(1,MAX_X+1,256))
    xs = sorted(xs)
    shots = 64
    sim = km.Simulator(batch_size=shots, counting_operations=True)
    sim.use_same_registers_as(circuit)
    scratch = arr(km.q(q) for q in range(22,circuit.num_qubits))
    worst = 0.0
    OUT.mkdir(parents=True,exist_ok=True)
    for checked,x in enumerate(xs,1):
        sim.clear_for_shot()
        sim.write_across_shots('x', ((np.array([x]) >> np.arange(22)) & 1).astype(bool)[:,None].repeat(shots,axis=1))
        sim.do(circuit)
        expected = scalar_horner(x) % 1
        actual_x = sim.read_across_shots('x')
        assert np.all(actual_x == (((x >> np.arange(22)) & 1).astype(bool)[:,None]))
        assert not sim.read_across_shots(scratch).any(), ('dirty scratch',x)
        for shot in range(shots):
            actual = Fraction(sim.read_shot_phase(shot),2) % 1
            assert actual == expected, ('phase mismatch',x,shot,actual,expected)
        rel = abs(float(scalar_horner(x))-1/math.sqrt(x))*math.sqrt(x)
        assert rel < 2**-24
        worst = max(worst,rel)
        if checked % 64 == 0:
            (OUT/'state.json').write_text(json.dumps({'phase':'testing','inputs':checked,'shots':checked*shots}))
    stats = circuit.operation_statistics()
    # The existing exhaustive integer probe supplies this all-input minimum slack.
    slack = json.loads((ROOT/'results/invsqrt-piecewise-math.json').read_text())['minimum_phase_slack_turns']
    finer = stats.num_z_pow(min_clifford_hierarchy_level=4)
    eps = math.pi*slack/finer
    rotation_t = finer*(0.56*math.log2(1/eps)+5.3)
    sampled = sim.operation_statistics()
    total_shots = len(xs)*shots
    mean_toffoli = sampled.num_toffoli_like/total_shots
    mean_t_like = sampled.num_t_like/total_shots
    c = mean_toffoli + (mean_t_like+rotation_t)/4
    return {'phase':'complete','sampled_circuit_passed':True,'inputs':len(xs),'shots':len(xs)*shots,
            'phase_check':'exact integer Horner phase on every sampled shot',
            'scratch_clean':True,'input_restored':True,'sample_max_relative_error':worst,
            'qubits':circuit.num_qubits,'toffoli_like':stats.num_toffoli_like,
            't_like':stats.num_t_like,'non_clifford_rotations':finer,
            'sample_mean_toffoli_like':mean_toffoli,'sample_mean_t_like':mean_t_like,
            'rotation_epsilon':eps,'rotation_t_gates_paper_formula':rotation_t,
            'toffoli_equivalents_paper_formula':c,'qubits_times_cost':c*circuit.num_qubits,
            'reference_qubits_times_cost':260*7988,'fraction_of_reference':c*circuit.num_qubits/(260*7988),
            'resource_gate_passed':c*circuit.num_qubits <= 0.5*260*7988,
            'official_verification':'pending official verifier and 200540928 shots',
            'operation_counts':{key:int(value) for key,value in stats.items()}}


if __name__ == '__main__':
    OUT.mkdir(parents=True,exist_ok=True)
    t0=time.process_time()
    (OUT/'state.json').write_text(json.dumps({'phase':'generating'}))
    try:
        circuit=make_circuit()
        circuit.to_file(OUT/'invsqrt.kmx')
        print(circuit.operation_statistics(),flush=True)
        result=probe(circuit)
        result['cpu_seconds']=time.process_time()-t0
        (OUT/'state.json').write_text(json.dumps(result,indent=2))
        print(json.dumps(result,indent=2))
    except Exception as e:
        (OUT/'failure.json').write_text(json.dumps({'phase':'failed','error':repr(e)}))
        raise
