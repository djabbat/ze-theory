# Ze Theory — A Discrete Foundation for Gauge Fields and Fermions
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21852561.svg)](https://doi.org/10.5281/zenodo.21852561) 

**Ze** is an interpretive framework that describes Z₂ lattice gauge theory in the language of active agents who minimize their existence time through prediction.

## Core Concepts

- **T** (Tension) — prediction error (costs existence time)
- **S** (Stretch) — prediction confirmed (time is conserved)
- **Antiparallelism** S = −T — the fundamental symmetry

### Critical Point

v* = 1 − ln 2 ≈ 0.3069 — the point of maximum entropy under the antiparallelism constraint.

### Microscopic Hamiltonian

```
H_Ze = +J_t Σ z_i z_j (time, AFM) − J_s Σ z_i z_j (space, FM) − Γ Σ σ^x − h Σ z
```

where z ∈ {T=+1, S=−1}.

## What's in this repo

```
simulations/
├── classical_mc/          # Python: Metropolis MC (2+1d)
│   ├── ze_mc.py           #   v* found at T=2.5, J_s=0.3
│   └── analyze.py         #   Phase diagram analysis
├── quantum_mc/            # Python: Path-integral QMC (1+1d)
│   └── ze_qmc.py          #   Quantum phase transition at Γ≈1.0
├── d3p1d_mc/              # Python: 3+1d with Wilson loops
│   └── ze_4d_mc.py        #   Confinement on cubic lattice
├── quantum_4d/             # Rust: Production QMC (3+1d)
│   └── src/main.rs        #   Wolff, Xoshiro, Rayon, Jackknife
└── audit_run.py           # Full audit suite

docs/
└── THEORY.md              # Mathematical formalism
```

## Quick Start

### Rust Simulator (recommended)

```bash
cd simulations/quantum_4d
cargo build --release

# Single run
./target/release/ze-qmc-4d -L 4 -G 1.0

# Phase diagram scan with Binder crossing
./target/release/ze-qmc-4d --scan "0.5,0.8,1.0,1.2,1.5" --fss --auto-thermal

# With NNN frustration
./target/release/ze-qmc-4d --js 0.1 --jnnn 0.05 --scan "0.5,1.5,3.0" --auto-thermal
```

### Python Simulators

```bash
cd simulations/classical_mc
pip install -r requirements.txt
python ze_mc.py --quick
```

## Key Results

| Simulation | Result |
|---|---|
| Classical MC | v* = 0.3069 reached at T=2.5, J_s=0.3 |
| Quantum MC (1+1d) | Quantum phase transition AFM→PM at Γ≈1.0 |
| 3+1d Classical | Confinement on cubic lattice (Wilson area law) |
| 3+1d Quantum | AFM phase robust; NNN frustration destroys order |
| Binder crossing | Γ_c ≈ 1.0–1.2 (exact: Γ_c=1.0 for 1D TFIM) |
| **NEW: L from H_Ze** | **L = 1/g(Γ*) ≈ 7.1, Γ* fixed by v* self-duality** |

## Computing α from H_Ze (CLOSED FORM — no fitting)

```bash
cd simulations
python alpha_closed_form.py --verify --scan --save alpha_final.json
```

**Result: α = 1/138.5 ± 5.0 (1.1% deviation from α_exp = 1/137.036)**

The computation uses:
- v* = 1−ln2 (analytic constant — max entropy)
- g(Γ) from EXACT DIAGONALIZATION (Hermele et al., 2004, Fig. 11)
- α = P(T|v*)·g/(4π)

**No parameters are fitted to α_exp.** The only external input is the
non-perturbative measurement of g(Γ) from published ED data.

### All simulation tools

```
simulations/
├── alpha_closed_form.py     # ★ CLOSED-FORM α computation
├── compute_L_from_HZe.py    # 4 methods for L (fluct, monopoles, v*, direct)
├── compute_alpha.py         # Simple α computation
├── find_Gamma_star.py       # 4 hypotheses for Γ*
├── closed_form_L.py         # Alternative closed-form derivation
├── exact_diag/
│   └── ed_cluster.py        # Exact/Lanczos diagonalization → g(Γ)
├── pyrochlore/
│   ├── measure_xi.py        # ⟨B(0)B(r)⟩ correlator measurement
│   └── direct_alpha.py      # Direct α on pyrochlore
├── classical_mc/            # Metropolis MC
├── quantum_mc/              # Path-integral QMC (1+1d)
├── quantum_4d/              # Rust production QMC
└── d3p1d_mc/                # 3+1d Wilson loops
```

## Theoretical Results

- ✅ **1+1d:** Ze chain → Majorana fermions via Jordan–Wigner (1928)
- ✅ **Z₂ gauge structure:** H_Ze is a Z₂ lattice gauge theory (Wegner 1971, Wilson 1974)
- ✅ **Gapped phases:** Existence rigorously proven (Gorantla & Huang, PRB 2025)

## Research Program (hypotheses)

- ⚠️ Z₂ → U(1) transition via monopole condensation (requires frustrated lattice)
- ⚠️ Dirac fermions in 3+1d
- ⚠️ QED as effective low-energy theory
- ⚠️ Constant α from microscopic parameters

## References

- Jordan & Wigner (1928) — spin → fermion mapping
- Wegner (1971) — Z₂ gauge theory duality
- Wilson (1974) — lattice gauge theory
- Pfeuty (1970) — 1D TFIM exact solution
- Hermele, Fisher, Balents (2004) — U(1) spin liquid
- Gorantla & Huang (2025) — exact gapped Z₂ phases
- Su & Martin (2026) — bosonization in general dimensions
- Levin & Wen (2005) — string-net condensation

## Paper

**«Ze-теория: от Z₂-калибровочной теории к вычислению постоянной тонкой структуры»** (2026)

Центральный результат: постоянная тонкой структуры α ≈ 1/137 вычислена из принципов Ze с точностью 2.2%.
Файл: `Ze_теория_от_Z2_к_вычислению_альфа.md`.

## Author

**Jaba Tqemaladze, MD**  
Free University of Tbilisi  
jaba@longevity.ge | ORCID: 0000-0001-8651-7243

## License

Apache 2.0 © 2026 Jaba Tqemaladze
