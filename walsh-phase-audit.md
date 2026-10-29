# A fixed Walsh-phase alternative

The accumulator-free idea is to approximate `g(t) = (1+t)^(-1/2)` on
21 normalized mantissa bits by diagonal Walsh terms, and tensor it with
`2^(-k/2)` on a five-bit exponent. This avoids the integer product register
but prices each emitted fine rotation.

## Why the tail sum matters

The Laplace representation

```text
g(t) = 1/sqrt(pi) * integral_0^infinity u^(-1/2) exp(-u) exp(-u*t) du
```

has nonnegative Walsh coefficients. For each input bit of weight `w`,
`exp(-u*w*b)` has coefficients `(1+exp(-u*w))/2` and
`(1-exp(-u*w))/2`, both nonnegative. Products and the positive integral
preserve this property. Discarding a coefficient set therefore incurs a
maximum absolute error equal to its sum, attained at mantissa `t=0`.

The fixed audit assigns half the arithmetic tolerance to truncation:
`2^-25 / sqrt(2) = 2.1073424e-8`. It evaluates 2,097,152 grid points.
The selected coefficient sum reserves **2e-11** for binary64 arithmetic:
21 butterfly stages and fewer than 3,000 selected terms.

| Fixed audit | Result |
|---|---:|
| Kept mantissa terms | 2,395 |
| Omitted sum estimate | 2.1034727e-8 |
| Omitted sum plus arithmetic allowance | 2.1054727e-8 |
| Nonzero exponent terms | 32 |
| Straight tensor rotations | 76,640 |
| Rotation-only paper price | 641,933 Toffoli-equivalents |
| Price times just the 22 input qubits | 14,122,533 |
| CPU time / peak RSS | 0.229 s / 180 MiB |

The paper pricing uses the remaining worst-input phase slack and
`.56 log2(1/eps) + 5.3` T gates per fine rotation. Normalization, controls,
and cleanup would add resources.

**Decision: stop this straight tensor route.** The same-day public inverse-root
resource comparison is about 104,000 qubit-Toffoli-equivalents, against
14.1 million from rotations alone here. This audit applies to the specified
truncation and exponent encoding. Different factorizations have their own
error and circuit contracts.

```sh
# Python, NumPy; Linux: 512 MiB address-space and 30 CPU-second limits.
python invsqrt-walsh-bound.py
```

The script writes to `temp/invsqrt-walsh-probe/`; `PROBE_OUT` selects a
different scratch directory. Its assertions cover the retained-positive
coefficients and the truncation budget with the arithmetic allowance.
[Measured output](results/invsqrt-walsh-bound.json). This is a mathematical
resource audit; native circuit construction stops at the price gate.

References: [benchmark paper](https://arxiv.org/abs/2610.23832),
[public resource comparison, 29 October 2026](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753527).
