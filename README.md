# Direct inverse-root phasing

AI engineering research into smaller, classically verifiable quantum circuits.

## First result: integer arithmetic

A **32-interval, degree-3 polynomial** directly approximates inverse square root, replacing the polynomial-seed-plus-Newton structure in the [Google Quantum AI benchmark](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation). Coefficients have 28 fractional bits; Horner products use integer floor shifts. Normalization uses the leading bit and folds odd-exponent scaling into the coefficient table.

The probe evaluated **all 3,133,452 valid integer inputs**:

| Measurement | Result |
|---|---:|
| Maximum relative error | 1.03631e-8 |
| Published arithmetic tolerance | 2^-24 ≈ 5.96046e-8 |
| Maximum fraction of tolerance used | 0.173864 |
| Tolerance headroom | 5.75× |
| Worst input | 999371 |
| Minimum phase slack | 2.90180e-11 turns |

The arithmetic contract passes. The stricter exploratory 8×-headroom check failed; both results are recorded in [the measured output](results/invsqrt-piecewise-math.json).

```sh
# Python 3.10+, NumPy; Linux resource limits: 1 GiB address space, 30 CPU seconds.
python invsqrt-piecewise-math.py
```

The script independently compares scalar and vectorized integer Horner values in every block, checks the complete input range, and writes coefficients and progress under `temp/invsqrt-probe/`. The measured run used 0.219 CPU seconds and about 74 MiB peak RSS.

## Next experiment

Synthesize the reversible kickmix circuit, account for scratch cleanup and rotation precision, and compare measured qubits × Toffoli-equivalents against the reference **260 × 7,988**. Gate cost and official verification are the next measurements.

References: [benchmark paper](https://arxiv.org/abs/2610.23832), [kickmix](https://github.com/quantumlib/Qualtran/tree/main/kickmix). The arithmetic probe uses the public mathematical contract and its own generated inputs.

Apache-2.0.
