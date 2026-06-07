# Part XXII — Multidimensional Tensor Convolution and the Spectral Cutoff of Lattice Transitions

*Volume II of the Arithmetic Geodynamics programme on the 6N skeleton.*

Parts XX–XXI worked at a single prime frequency q on the torus `T¹`. This paper
takes the product over **all** prime frequencies — the space `∏_{q>3} T¹_q` — and
studies the global twin–twin correlation

> R(j) = ∏_{q>3} R_q(j).

**Main result: R(j) converges absolutely to a strictly positive constant**, with a
sharp mechanism:

- **Resonance is finite.** The constructive states deviate by O(1/q) but occur at
  only finitely many primes:
  - State A `R_q = q/(q−2)` requires `q | j`;
  - State B `R_q = q(q−3)/(q−2)²` requires `q | 3j±1`.
  All such primes are ≤ 3j+1.
- **The tail is locked.** Every prime beyond is State C,
  `R_q = q(q−4)/(q−2)² = 1 − 4/(q−2)² = O(1/q²)`, summable — so the product
  converges (a quasi-periodic diffraction pattern dominated by small primes).

### Spectral cutoff k\*(ε)

Full product (j=1, all State C): `R(1) = 𝔖_quad/𝔖_twin² = 0.3969`.

| ε | k\* | primes retained |
|---|---|---|
| 10⁻¹ | 13 | 4 |
| 10⁻² | 73 | 19 |
| 10⁻³ | 557 | 100 |
| 10⁻⁴ | 4229 | 577 |

### The Spatial Transition Operator

`𝓜_{A,B}(j) = ρ_B · ∏_q R_{q,A,B}(j)` — the Hardy–Littlewood joint singular series
re-expressed as a conditional density (state B at N+j given state A at N). It is a
quasi-periodic, long-range, small-prime-dominated correlation — a multi-frequency
diffraction grating — exact at integer lags and equal there to the discrete sieve
ratio. It is a compact bookkeeping of inter-centre correlation, **not** a dynamical
law.

## Layout

```
.
├── paper/    Chen_6N_Paper22.{tex,pdf} + figure
├── figures/  fig_spectral_cutoff.{pdf,png}
├── data/     kstar_table.csv · R1_product.csv
├── code/
│   ├── fig_spectral_cutoff_make.py   # resonance-vs-tail + k*(ε) figure
│   └── verify_spectral_cutoff.py     # resonance-finiteness checks; writes data/
├── CITATION.cff · .zenodo.json · LICENSE (MIT)
```

## Reproducing

```bash
pip install numpy matplotlib
python code/verify_spectral_cutoff.py   # checks resonance finiteness; (re)writes data/, prints k*(ε)
python code/fig_spectral_cutoff_make.py # regenerates the figure
```

Expected: resonance-finiteness checks pass; `R(1)=0.396881`; the k\*(ε) table above.

## Scope

Closed-form throughout. At integer lags the infinite product equals the discrete
singular-series ratio exactly; no new prime data, no infinitude claim. Continues
Part XXI (doi:10.5281/zenodo.20584909).

## License

MIT — see `LICENSE`.
