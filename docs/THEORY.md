# Ze Theory — Mathematical Formalism

## 1. Axioms

1. **Binarity:** Fundamental reality is a binary variable z ∈ {T=+1, S=−1}
2. **Antiparallelism:** S = −T (T and S are opposite)
3. **Locality:** Spatially close systems have correlated T/S states

## 2. Ze Parameters

| Parameter | Formula | Meaning |
|-----------|---------|---------|
| **v** (Ze velocity) | (N_T − N_S) / N | Balance of T/S; v* = 1−ln 2 ≈ 0.3069 |
| **τ** (Ze complexity) | H(stream) / log₂(N) | Normalized entropy |
| **Z** (Ze index) | N_T / N | Fraction of T-events |

## 3. Hamiltonian

```
H_Ze = -J_t Σ_{x,t} z_{x,t} z_{x,t+1}           [time: antiparallelism]
       -J_s Σ_{⟨x,y⟩,t} z_{x,t} z_{y,t}          [space: locality]
       -h Σ_{x,t} z_{x,t}                         [field: T/S asymmetry]
```

where J_t > 0 (antiferromagnetic), J_s > 0 (ferromagnetic).

## 4. Z₂ Gauge Structure

H_Ze possesses local Z₂ symmetry: flipping z at one site and all incident spatial links leaves the Hamiltonian invariant. This makes H_Ze a Z₂ lattice gauge theory (Wegner 1971).

**Gauss law operator:**
```
G_x = Π_μ σ^x_{x,μ}
```
Physical states satisfy G_x|Ψ⟩ = +|Ψ⟩.

## 5. Phase Diagram

Z₂ gauge theory in 3+1d has two phases:
- **Confinement** (strong coupling): gapped, electric charges confined
- **Deconfinement** (weak coupling): topological order, degenerate ground states

At the frustration-free point (Gorantla & Huang 2025): exactly 9-fold degenerate ground states with proven gap.

## 6. Path to QED

1. ✅ Z₂ gauge theory from Ze axioms (Wegner 1971)
2. ✅ 1+1d Majorana fermions via Jordan–Wigner (1928)
3. ⚠️ Z₂ → U(1) transition via monopole condensation
4. ⚠️ Dirac fermions in 3+1d
5. ⚠️ QED as effective low-energy theory

## 7. Key Papers

- Wegner, F.J. (1971). *J. Math. Phys. 12, 2259*
- Wilson, K.G. (1974). *Phys. Rev. D 10, 2445*
- Jordan, P. & Wigner, E. (1928). *Z. Phys. 47, 631*
- Gorantla, P. & Huang, T.-C. (2025). *Phys. Rev. B 111, 245110*
- Su, L. & Martin, I. (2026). *SciPost Phys. 20, 180*
- Levin, M.A. & Wen, X.-G. (2005). *Phys. Rev. B 71, 045110*
- Hermele, M. et al. (2004). *Phys. Rev. B 69, 064404*
