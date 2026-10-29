#!/usr/bin/env python3
"""One fixed Walsh phase audit; no circuit/native/parameter scan.

f(t)=1/sqrt(1+t) has nonnegative Walsh coefficients: its Laplace
representation is a positive integral of products of exp(-u*b_i*2^-i),
and each binary factor has two nonnegative Walsh coefficients. Therefore
a discarded coefficient set has worst error equal to its sum at t=0.
The audit covers a sufficient uniform-mantissa error rule, not all possible
phase approximations. A different exponent encoding is a separate design.
"""
import json, math, time, os, resource
from pathlib import Path
os.environ['OPENBLAS_NUM_THREADS'] = '1'
resource.setrlimit(resource.RLIMIT_CPU, (30, 30))
resource.setrlimit(resource.RLIMIT_AS, (512 * 1024**2, 512 * 1024**2))
import numpy as np
out=Path(os.environ.get('PROBE_OUT', 'temp/invsqrt-walsh-probe'))
out.mkdir(parents=True, exist_ok=True)
t=time.process_time()
n=21
N=1<<n
v=1/np.sqrt(1+np.arange(N,dtype=np.float64)/N)
h=1
while h<N:
    p=v.reshape(-1,2*h)
    a=p[:,:h].copy()
    b=p[:,h:].copy()
    p[:,:h]=a+b
    p[:,h:]=a-b
    h*=2
v/=N
ranked=np.sort(v)[::-1]
prefix=np.cumsum(ranked,dtype=np.float64)
# Half of relative 2^-24 assigned to truncation, at min normalized f=1/sqrt(2).
B=2**-25/math.sqrt(2)
# The unnormalized FWHT has 21 addition/subtraction stages.  The selected
# coefficient sum has at most K*gamma_21 plus gamma_K absolute error here.
# 2e-11 covers that bound for K < 3000 and bounded IEEE binary64 sqrt error.
roundoff=2e-11
keep=int(np.searchsorted(prefix, 1-B+roundoff, side='left'))+1
assert 0 < keep < 3000 and np.all(ranked[:keep] > 0)
omitted=1-float(prefix[keep-1])
assert 0 <= omitted and omitted+roundoff <= B
assert 1-float(prefix[keep-2])+roundoff > B
# Direct binary exponent factor, h(k)=2^(-k/2), including unused k=22..31.
e=2**(-np.arange(32,dtype=np.float64)/2)
h=1
while h<32:
    p=e.reshape(-1,2*h)
    a=p[:,:h].copy(); b=p[:,h:].copy()
    p[:,:h]=a+b; p[:,h:]=a-b; h*=2
e/=32
# Relevant uniform phase slack if above arithmetic bound is met.
slack=2**-25/math.sqrt(3133452)
rot=keep*32
teq_per=(.56*math.log2(rot/(math.pi*slack))+5.3)/4
report={
 'basis':'Walsh phase expansion on 21-bit normalized mantissa',
 'mantissa_grid_points':N,
 'truncation_l1_budget':B,
 'kept_terms_for_uniform_l1_rule':keep,
 'omitted_l1_sum_estimate':omitted,
 'omitted_l1_with_roundoff_allowance':omitted+roundoff,
 'largest_omitted_coefficient_estimate':float(ranked[keep]),
 'previous_omitted_sum_estimate':1-float(prefix[keep-2]),
 'total_coefficient_sum_estimate':float(prefix[-1]),
 'exponent_basis_terms':int(np.count_nonzero(e)),
 'exponent_coefficient_abs_range':[float(np.min(np.abs(e))),float(np.max(np.abs(e)))],
 'straight_tensor_rotations':rot,
 'half_tolerance_phase_slack_at_max_x':slack,
 'straight_tensor_rotation_only_cost':rot*teq_per,
 'straight_tensor_rotation_price_times_22_input_qubits':22*rot*teq_per,
 'floating_audit_roundoff_allowance':roundoff,
 'scope':'fixed top-coefficient truncation and 32-term exponent tensor; different encodings priced separately',
 'native_circuit_built':False,
 'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
 'cpu_seconds':time.process_time()-t,
}
(out/'invsqrt-walsh-bound.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
