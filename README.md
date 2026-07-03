# Ze Theory — A Discrete Foundation for Gauge Fields and Fermions

**Ze** is a theory that derives quantum electrodynamics from a single postulate: reality is a binary stream of T/S events.

- **T** (Tension) — prediction error (costs proper time)
- **S** (Stretch) — prediction confirmed (time stands still)
- **Antiparallelism** S = −T — the fundamental symmetry

### Critical Point

v* = 1 − ln 2 ≈ 0.3069 — the point of maximum entropy with antiparallelism constraint.

### Microscopic Hamiltonian

```
H_Ze = -J_t Σ z_i z_j (time) - J_s Σ z_i z_j (space) - h Σ z_i
```

where z ∈ {T=+1, S=−1}.

### What's in this repo

```
simulations/classical_mc/
├── ze_mc.py          # Metropolis Monte Carlo (Numba-optimized)
├── analyze.py        # Phase diagram analysis
├── requirements.txt
└── README.md

docs/
└── THEORY.md         # Mathematical formalism
```

### Quick Start

```bash
cd simulations/classical_mc
pip install -r requirements.txt
python ze_mc.py --quick
```

### Theoretical Results (proven)

- ✅ **1+1d:** Ze chain → Majorana fermions via Jordan–Wigner (1928)
- ✅ **Z₂ gauge structure:** H_Ze is a Z₂ lattice gauge theory (Wegner 1971, Wilson 1974)
- ✅ **Gapped phases:** Existence rigorously proven (Gorantla & Huang, PRB 2025)

### Research Program (hypotheses)

- ⚠️ Z₂ → U(1) transition via monopole condensation
- ⚠️ Dirac fermions in 3+1d
- ⚠️ QED as the effective low-energy theory

### References

- Jordan & Wigner (1928) — spin → fermion mapping
- Wegner (1971) — Z₂ gauge theory duality
- Wilson (1974) — lattice gauge theory
- Su & Martin (2026) — bosonization in general dimensions
- Gorantla & Huang (2025) — exact gapped Z₂ phases
- Levin & Wen (2005) — string-net condensation

### Author

**Jaba Tqemaladze, MD**  
Free University of Tbilisi  
jaba@longevity.ge | ORCID: 0000-0001-8651-7243

### License

Apache 2.0 © 2026 Jaba Tqemaladze
