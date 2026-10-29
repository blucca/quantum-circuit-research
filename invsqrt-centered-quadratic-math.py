"""Frozen centered quadratic witness, 128 cells / 29 coefficient frac bits.
Only mathematical error is measured here; no competitive circuit is claimed.
The exact signed local r is the 14-bit low word with its top bit flipped.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json, math, resource, time
from pathlib import Path
resource.setrlimit(resource.RLIMIT_CPU,(20,20))
resource.setrlimit(resource.RLIMIT_AS,(512*1024**2,512*1024**2))
import numpy as np
OUT=Path(os.environ.get("PROBE_OUT", "temp/invsqrt-centered-quadratic-math-probe"))
OUT.mkdir(parents=True, exist_ok=True)
N,CELLS,FRAC,LOW,CHUNK=3_133_452,128,29,14,1<<18
center=1+(np.arange(CELLS,dtype=np.float64)+.5)/CELLS
# t=r/2^13 in [-1,1); mantissa displacement is t/256.
f=1/np.sqrt(center)
b=-.5*center**-1.5/256
c=.375*center**-2.5/256**2
coef=np.rint(np.stack([np.stack([f,b,c],axis=-1)/math.sqrt(2**o)
                       for o in (0,1)])*2**FRAC).astype(np.int64)
assert np.all((0<=coef[...,0]) & (coef[...,0]<2**29))
assert np.all((-2**20<=coef[...,1]) & (coef[...,1]<0))
assert np.all((0<=coef[...,2]) & (coef[...,2]<2**12))
(OUT/'coefficients.json').write_text(json.dumps(coef.tolist()))
def scalar(x):
 k=int(x).bit_length()-1;m=int(x)<<(21-k)
 s=(m>>LOW)&127;r=(m&((1<<LOW)-1))-(1<<13)
 a,b,c=map(int,coef[k&1,s]);return (a<<26)+(b*r<<13)+c*r*r,55+k//2
start=time.process_time();maximum=0.;slack=float('inf');state={}
for first in range(1,N+1,CHUNK):
 x=np.arange(first,min(first+CHUNK,N+1),dtype=np.int64)
 k=np.floor(np.log2(x)).astype(np.int64)
 assert np.all((1<<k)<=x) and np.all(x<(1<<(k+1)))
 m=np.left_shift(x,21-k);s=(m>>LOW)&127;r=(m&((1<<LOW)-1))-(1<<13)
 z=coef[k&1,s]
 numerator=(z[:,0]<<26)+((z[:,1]*r)<<13)+z[:,2]*r*r
 assert np.all((0<numerator)&(numerator<2**56))
 for i in (0,len(x)//2,len(x)-1):assert scalar(int(x[i]))==(int(numerator[i]),int(55+k[i]//2))
 estimate=np.ldexp(numerator.astype(np.float64),-55-k//2);exact=1/np.sqrt(x.astype(np.float64))
 error=np.abs(estimate-exact);relative=error/exact;j=int(relative.argmax())
 if relative[j]>maximum:maximum=float(relative[j]);worst=int(x[j])
 slack=min(slack,float((exact*2**-24-error).min()))
 state={'phase':'math-running','checked':int(x[-1]),'maximum_relative_error':maximum,'worst_x':worst,'minimum_phase_slack_turns':slack}
 (OUT/'state.json').write_text(json.dumps(state,indent=2))
state.update(phase='math-complete',math_passed=maximum<2**-24 and slack>0,
 tolerance_ratio=maximum*2**24,cpu_seconds=time.process_time()-start,
 coefficient_bits=[29,21,12],local_signed_bits=14,exact_numerator_bits=56,
 competitive_circuit='not built; full lookup/square/cleanup must pass current 104k gate')
(OUT/'state.json').write_text(json.dumps(state,indent=2));print(json.dumps(state,indent=2))
assert state['math_passed'],'Frozen quadratic fails exact all-input contract'
