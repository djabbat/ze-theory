# Ze Theory — Mathematical Formalism

## 1. Axioms

1. **Binarity:** Fundamental reality is a binary variable z ∈ {T=+1, S=−1}
2. **Antiparallelism:** S = −T (T and S are opposites)
3. **Locality:** Spatially close systems have correlated T/S states

## 2. Ze Parameters

| Parameter | Formula | Meaning |
|-----------|---------|---------|
| **v** (Ze velocity) | (N_T − N_S) / N | Balance of T/S events |
| **v*** (critical) | 1 − ln 2 ≈ 0.3069 | Maximum entropy with antiparallelism |
| **τ** (Ze complexity) | H(stream) / log₂(N) | Normalized entropy |
| **τ_Ze** (existence time) | η · N_T | Resource depleted by T-events |

## 3. Hamiltonian

```
H_Ze = +J_t Σ z_i z_j (time: antiferromagnetic, antiparallelism)
       −J_s Σ z_i z_j (space: ferromagnetic, locality)
       −Γ Σ σ^x        (quantum fluctuations)
       −h Σ z          (external field)
```

Quantum Monte Carlo via Suzuki-Trotter decomposition maps this to a classical model in
(d+1) dimensions with effective Trotter coupling K_τ = −½ ln tanh(βΓ/M).

## 4. Z₂ Gauge Structure

H_Ze possesses local Z₂ symmetry: flipping z at one site and all incident spatial links
leaves the Hamiltonian invariant. This makes H_Ze a Z₂ lattice gauge theory (Wegner 1971).

**Gauss law operator:** G_x = σ^x_x Π_μ σ^x_{x,μ}.
Physical states satisfy G_x|Ψ⟩ = +|Ψ⟩.

## 5. Phase Diagram

Z₂ gauge theory in 3+1d has three phases (Wilson 1974, Fradkin 2013):

| Phase | Parameters | Properties | Ze Interpretation |
|-------|-----------|------------|-------------------|
| **Confinement** | Γ ≪ J | Gapped, area law for Wilson loops | T/S events strongly correlated; AFM order |
| **Deconfinement** | Γ ≫ J | Topological order, perimeter law | T/S events independent; gauge field liberated |
| **Higgs** | h ≫ J | Fully symmetric | External field fixes T/S bases |

On a cubic lattice, only confinement and Higgs phases are observed.
U(1) spin liquid (deconfinement) requires geometrical frustration (pyrochlore lattice).

## 6. Path to QED

1. ✅ Z₂ gauge theory from Ze axioms (Wegner 1971)
2. ✅ 1+1d Majorana fermions via Jordan–Wigner (1928)
3. ⚠️ Z₂ → U(1) transition via monopole condensation (requires frustration)
4. ⚠️ Dirac fermions in 3+1d
5. ⚠️ QED as effective low-energy theory

## 7. Numerical Methods

- **Classical MC:** Metropolis algorithm on Lx×Ly×Lt lattice
- **Quantum MC:** Path-integral with Wolff cluster updates
- **Observables:** energy, magnetization, staggered magnetization, Binder cumulant, Wilson loops
- **Error analysis:** Jackknife resampling, integrated autocorrelation time τ_int
- **Performance:** i8 compressed storage (8× memory savings), Rayon parallelism, Xoshiro256 RNG

## 8. Connection to FEP

The Free Energy Principle (Friston 2010) states that systems minimize variational free energy
(surprisal). In the binary Ze framework, F(v) = −ln P(S) − H(v) is minimized at v = v*.
Ze provides a microscopic, binary realization of the FEP variational principle.

## 9. Key Papers

- Jordan, P. & Wigner, E. (1928). *Z. Phys.* 47, 631 — spin→fermion mapping
- Wegner, F.J. (1971). *J. Math. Phys.* 12, 2259 — Z₂ gauge duality
- Wilson, K.G. (1974). *Phys. Rev. D* 10, 2445 — lattice gauge theory
- Pfeuty, P. (1970). *Ann. Phys.* 57, 79 — 1D TFIM exact solution
- Hermele, M. et al. (2004). *Phys. Rev. B* 69, 064404 — pyrochlore photons
- Gorantla, P. & Huang, T.-C. (2025). *Phys. Rev. B* 111, 245110 — gapped Z₂ phases
- Su, L. & Martin, I. (2026). *SciPost Phys.* 20, 180 — bosonization in general dimensions
- Levin, M.A. & Wen, X.-G. (2005). *Phys. Rev. B* 71, 045110 — string-net condensation
- Friston, K. (2010). *Nat. Rev. Neurosci.* 11, 127 — free energy principle
- Fields, C. et al. (2022). *Prog. Biophys. Mol. Biol.* 173, 36 — FEP for quantum systems

## 10. Proper Time Synchrony Principle (PTSP) — NEW 2026-07-28

### 10.1 Theorem

If N subsystems are created with identical proper time τᵢ(t₀) = τ₀ and identical proper time flow rate dτᵢ/dt|t₀ = ω₀, then:

1. **While** τᵢ(t) = τⱼ(t) for all i,j: **dS/dt = 0** (entropy conserved)
2. **Once** ∃i,j: τᵢ(t) ≠ τⱼ(t): **dS/dt = κ · Var(τ) + O(Var²)**

where Var(τ) = (1/N) Σᵢ (τᵢ − ⟨τ⟩)², κ > 0.

### 10.2 Proof (Ze framework)

Proper time for subsystem i: τᵢ = η · N_T⁽ⁱ⁾, where N_T⁽ⁱ⁾ is the count of T-events (prediction errors).

When N_T⁽ⁱ⁾ = N_T⁽ʲ⁾ for all i,j, all conditional distributions p(zᵢ | N_T) are identical → joint entropy factorizes → dS/dt = 0.

When N_T values diverge, p(zᵢ | N_T⁽ⁱ⁾) ≠ p(zⱼ | N_T⁽ʲ⁾) → additional microstates → dS/dt > 0, proportional to variance of τ.

### 10.3 Three Applications

| Domain | Meaning of τ | Consequence of Var(τ) > 0 |
|--------|-------------|---------------------------|
| **GR/Thermodynamics** | Proper time along worldlines | Gravity creates entropy; 2nd law = geometric effect |
| **Ze-Hierarchy** | Bot battery voltage V_cap | Hierarchy emerges (HI ∝ Var(V_cap)) |
| **MCARA/CEDAR** | Cellular counters (centriole, epigenetic, etc.) | Tissue aging; CEDAR resets Var(τ) → 0 for rejuvenation |

### 10.4 Key Equation

> dS(Σ)/dt = ∫ dμ(x) ∫ dμ(y) K(x,y) · |τ(x) − τ(y)|²

Full document: `Ze_Model/PROPER_TIME_ENTROPY.md`
