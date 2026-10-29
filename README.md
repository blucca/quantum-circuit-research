# Direct inverse-root phasing

AI engineering research into smaller, classically verifiable quantum circuits.

## Latest result: one-product biased linear phase

A **1024-cell linear interpolant** uses half the local chord-error bound as a
bias. Its joint QROM stores a 29-bit intercept and a 19-bit local cell change.
One exact signed product produces the positive 40-bit numerator
`(intercept << 11) + local_change * delta`; normalization and phase scaling
are reversed, and every scratch register is cleared.

| Measurement | Result |
|---|---:|
| Integer inputs checked | All 3,133,452 |
| Maximum relative error | 4.60454e-8 |
| Published arithmetic tolerance | 2^-24 ≈ 5.96046e-8 |
| Native circuit checks | 482 inputs × 64 HMR shots = 30,848 passing shots |
| Qubits | 150 |
| Mean Toffoli-like gates | 4,226 |
| T gates / finer-than-T rotations | 1 / 47 |
| Toffoli-equivalents, paper synthesis formula | 4,555.28 |
| Qubits × cost | 683,291 |
| Fraction of the paper reference | **0.329× — a 3.04× improvement** |

The sampled circuit restores every input, clears every scratch wire, and
matches the integer numerator's phase exactly. Arithmetic took 0.185 CPU
seconds; circuit generation and sampled validation took 0.766 CPU seconds.
The [arithmetic result](results/invsqrt-linear-math.json),
[joint coefficients](results/invsqrt-linear-coefficients.json), and
[native resource report](results/invsqrt-linear-circuit-kickmix.json) are available.
Validation covers exhaustive integer arithmetic and sampled native execution;
the official complete 200,540,928-shot verifier is a further stage.

```sh
python invsqrt-linear-math.py
# Install public kickmix as described below, then:
PROBE_OUT=temp/invsqrt-linear-circuit-probe python invsqrt-linear-kickmix.py
```

The frozen circuit uses 150 qubits, with only three retained clean wires during
phase scaling. The ten zero wires are released before multiplication cleanup.
The linear circuit reuses the leading-bit helper from `invsqrt-kickmix.py`.

## Alternative checked: accumulator-free Walsh phases

A [fixed Walsh audit](walsh-phase-audit.md) retains 2,395 normalized-mantissa
terms under a uniform error rule. Its straight exponent tensor needs 76,640
fine rotations; rotation-only pricing already exceeds the current public
comparison by about 136×. The audit and a small runnable script are public.

**Current research status:** arithmetic and sampled-circuit witnesses are
complete. Active optimization and formal entry are paused: the latest linear
circuit improves the paper reference by 3.04×, while the public same-task
comparison reports about 20×.

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

## Measured reversible circuit

[`invsqrt-kickmix.py`](invsqrt-kickmix.py) implements the same integer arithmetic
using public kickmix: leading-bit normalization, coefficient QROM, signed
multiply-adds, exponent-dependent phasing, and exact scratch cleanup. It retains
product history for a transparent first implementation.

**565 boundary and random inputs × 64 measurement shots = 36,160 passing shots.**
Every sampled input is restored, every scratch wire is zero, and every phase
matches the integer Horner result exactly. The probe took 0.562 CPU seconds.

| Resource | Measured value |
|---|---:|
| Qubits | 284 |
| Mean Toffoli-like gates | 9,252 |
| T gates | 1 |
| Finer-than-T rotations | 35 |
| Toffoli-equivalents, paper synthesis formula | 9,487.19 |
| Qubits × cost | 2,694,361 |
| Fraction of the paper reference | 1.297× |

[Measured output](results/invsqrt-circuit-kickmix.json) includes operation counts
and rotation pricing. The existing full-input arithmetic probe supplies the
minimum phase slack used to price synthesis. The complete 200,540,928-shot
competition verifier remains a separate validation stage.

```sh
# Install kickmix from the public Qualtran kickmix directory, then:
PROBE_OUT=temp/invsqrt-circuit-probe python invsqrt-kickmix.py
```

The research gate targeted half the reference resource cost. **This first
cubic circuit failed that resource gate; active optimization of this recipe is
stopped.** The integer approximation and reproducible circuit remain available
as research artifacts. Competition submission: none.

References: [benchmark paper](https://arxiv.org/abs/2610.23832),
[kickmix](https://github.com/quantumlib/Qualtran/tree/main/kickmix), public source
commit `638777fac4379bc62e2a1e4a52cf1a1753895d48`.

Apache-2.0.
