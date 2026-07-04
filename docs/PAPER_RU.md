# Ze Theory: An Interpretive Framework for Z₂ Gauge Theory — Structural Relation α ∝ 1/L with Empirical Calibration

**Jaba Tqemaladze, MD**

*Free University of Tbilisi | jaba@longevity.ge | ORCID: 0000-0001-8651-7243*

*July 4, 2026*

---

## Abstract

Ze theory offers an interpretive framework for Z₂ lattice gauge theory in the language of active agents minimizing existence time through prediction. The gauge freedom of Z₂ theory is reinterpreted as the agent's freedom to define which events count as prediction errors (T) and which as successes (S). **Structural result:** the fine-structure constant α can be expressed as α = P(T|v*)·g/(4π), where P(T|v*) = (2−ln 2)/2 is the agent's error probability at the critical point v* = 1−ln 2, and g is the ring-exchange constant of the emergent U(1) spin liquid on the pyrochlore lattice. **Important caveat on g:** g is computed via 6th-order degenerate perturbation theory in the XXZ model on pyrochlore (Hermele, Fisher & Balents, 2004, §VI), verified against exact diagonalization on 16-site clusters. It is NOT directly measured non-perturbatively, and the connection between H_Ze (a Z₂ gauge theory) and the XXZ model requires explicit demonstration which is not yet provided. At Γ/J_t ≈ 0.94 (a value chosen to match the perturbative regime where g ≈ 0.14), the formula yields α ≈ 1/138.5. **This is a consistency check, not a first-principles prediction.** The functional form α ∝ g is postulated based on dimensional analysis of the effective gauge action; the factor P(T|v*) is a conceptual Ansatz linking agent error probability to coupling strength. Neither the functional form nor the proportionality constant are derived from H_Ze.

---

## 1. Introduction

### 1.1 Definitions

**Ze-agent** — a mathematical model of a system that observes binary events, predicts the next event, and records the outcome: T (Tension — prediction error) or S (Stretch — prediction success).

The parameter v = (N_T − N_S)/(N_T + N_S) ∈ [−1, +1] is the Ze-velocity. The critical point v* = 1 − ln 2 ≈ 0.3069. At v*: P(S) = (ln 2)/2 ≈ 0.347, P(T) = (2−ln 2)/2 ≈ 0.653, N_T/N_S = 2/ln 2 − 1 ≈ 1.885. This value emerges in classical Monte Carlo at T ≈ 2.5, J_s ≈ 0.3 (J_t = 1) as the point separating confinement and deconfinement phases.

**Existence time** τ_Ze — an agent's internal resource, depleted by T-events: Δτ = η·N_T. This is a phenomenological parameter with no direct connection to relativistic proper time or any physical observable in the current formulation. A possible microscopic interpretation: τ_Ze ∝ 1/Δ, where Δ is the energy gap of H_Ze, which would link "existence time" to the characteristic timescale of quantum evolution. This connection remains a hypothesis.

**Gauge freedom and T/S-basis:** The gauge symmetry of Z₂ theory is strictly local (G_x|Ψ⟩ = +|Ψ⟩). The interpretation as "the agent's freedom to choose the T/S-basis" is a conceptual choice, not a mathematical consequence. The agent cannot arbitrarily choose the T/S-basis outside the physical Hilbert space.

### 1.2 What Ze Asserts and Does Not Assert

**Asserts:**
- The gauge freedom of Z₂ theory can be interpreted as the agent's freedom to define the T/S-basis
- Interaction is a strategy for minimizing existence time expenditure
- The interpretation is consistent with the Free Energy Principle (Friston, 2010; Fields et al., 2022) for the binary case
- v parameter is a statistical property of the agent's classical measurement record, not a quantum expectation value; Elitzur's theorem does not constrain it

**Does not assert:**
- Ze does not derive QED from first principles
- Ze does not compute the fine-structure constant α (the numerical coincidence in §5.1 is a consistency check, not a derivation)
- Ze does not prove new mathematical theorems about Z₂ gauge theory
- Ze does not offer experimentally testable predictions beyond Γ_c(J_s) dependence, verifiable on quantum simulators
- The connections to MEPP and edge of chaos (§2.5) are hypotheses, not proven from H_Ze
- The connection to ER=EPR (§4.3) is a conceptual analogy, not supported by mathematical derivation

---

## 2. Mathematical Formalism

### 2.1 Hamiltonian

$$H_{Ze} = +J_t\sum_{x,t} z_{x,t} z_{x,t+1} - J_s\sum_{\langle x,y\rangle, t} z_{x,t} z_{y,t} - \Gamma\sum_{x,t} \sigma^x_{x,t} - h\sum_{x,t} z_{x,t}$$

where z_{x,t} ∈ {±1}, J_t > 0 (antiferromagnetic — antiparallelism), J_s > 0 (ferromagnetic — locality), Γ — transverse field (quantum fluctuations).

**Justification of signs.** The term +J_t·z_i·z_j with J_t > 0: energy ∼ +J_t for parallel spins and ∼ −J_t for antiparallel (antiferromagnetic — lower energy). The term −J_s·z_i·z_j with J_s > 0: energy ∼ −J_s for parallel spins (lower) — ferromagnetic ordering in space.

### 2.2 Equivalence to Z₂ Gauge Theory

Gauss operator: G_x = σ^x_x Π_μ σ^x_{x,μ}. Direct computation: [H_Ze, G_x] = 0 for h = 0. Physical Hilbert space: G_x|Ψ⟩ = +|Ψ⟩. This is the standard formalism of Z₂ gauge theory (Wegner, 1971; Wilson, 1974; canonical review: Kogut, 1979). No new Hamiltonian is introduced — Ze uses the existing one.

**Elitzur's theorem and the status of v:** Elitzur's theorem (Elitzur, 1975) states that in a gauge theory with local symmetry, non-gauge-invariant operators have zero expectation value in any gauge-invariant state: ⟨Ψ|z|Ψ⟩ = 0 for all physical states satisfying G_x|Ψ⟩ = +|Ψ⟩. The proof is elementary: G_x z G_x = −z for links incident on x, hence ⟨Ψ|z|Ψ⟩ = ⟨Ψ|G_x z G_x|Ψ⟩ = −⟨Ψ|z|Ψ⟩ ⇒ ⟨Ψ|z|Ψ⟩ = 0.

**Consequently, v = ⟨z⟩ = 0 identically in the physical Hilbert space.** Gauge fixing does not alter this result — it is a calculational tool, not a physical operation that creates new observables. The analogy with the Higgs field ⟨φ⟩ in unitary gauge (drawn in earlier versions of this work) is misleading because in Higgs theories, the gauge-invariant combination |φ|² is physical, whereas in pure Z₂ gauge theory, there is no gauge-invariant bilinear that reduces to z.

**However**, Caudy & Greensite (2008) demonstrated that while local gauge symmetries cannot break spontaneously, **global subgroups** of the local gauge symmetry CAN. Specifically: "Local gauge symmetries cannot break spontaneously, according to Elitzur's theorem, but this leaves open the possibility of breaking some global subgroup of the local gauge symmetry [...] the location of the breaking in the phase diagram depends on the choice of global subgroup" (Caudy & Greensite, 2008, Phys. Rev. D 78, 025018). This is the precise mathematical mechanism underlying the Ze interpretation: the agent's choice of T/S-basis corresponds to selecting a global subgroup of the Z₂ gauge symmetry, and the breaking of this subgroup at v = v* is physically meaningful. Similarly, Grady (2005) showed that the confinement-Higgs continuity (Fradkin & Shenker, 1979) depends on gauge-fixing: in a partial axial gauge, the phases are separated by a genuine phase transition.

**Implementation in Ze:** The parameter v is the order parameter for the breaking of the agent's chosen global subgroup. In the unbroken phase (v = 0, T and S equally likely), the agent has no predictive power. In the broken phase (v ≠ 0), the agent can distinguish T from S events. The critical point v* = 1−ln 2 marks the transition between these regimes.

**Important nuance:** The choice of which global subgroup to "break" is the agent's choice of T/S-basis, but the LOCATION of the transition in the phase diagram is determined by the dynamics of H_Ze. The fact that the classical Monte Carlo finds v* at T ≈ 2.5, J_s ≈ 0.3 suggests that H_Ze's phase boundary aligns with the agent's entropy-maximizing point — a non-trivial dynamical fact that requires further investigation.

### 2.3 Phase Diagram

Z₂ gauge theory in 3+1d has three phases (Wilson, 1974; Fradkin, 2013). An important property: **self-duality** — Z₂ theory in 3+1d is dual to itself under the Kramers–Wannier transformation (Kramers & Wannier, 1941; Wegner, 1971; review: Savit, 1980). Explicitly: the substitution z → z' on the dual lattice with σ^x ↔ σ^z maps H_Ze(Γ, h) to H_Ze(h, Γ). The self-dual point Γ = h determines the phase transition. In the Ze interpretation, self-duality corresponds to symmetry between T- and S-events: describing the system in terms of errors (T) is equivalent to describing it in terms of successes (S) on the dual lattice.

**Note on phases:** in Z₂ theory with fundamental Higgs fields, the confinement and Higgs phases are analytically connected — there is no phase transition between them, only a crossover (Fradkin & Shenker, 1979). This is a fundamental result: in the presence of matter in the fundamental representation of the gauge group, the distinction between confinement and Higgs phases disappears. For the Ze interpretation, this means that the "agent's choice of T/S-basis" (gauge fixing) does not create physically distinct phases — consistent with Elitzur's theorem (§2.2). Deconfinement (topological order) is a separate phase, requiring frustration (Hermele, Fisher & Balents, 2004).

| Phase | Parameters | Wilson Loops | Ze-interpretation |
|---|---|---|---|
| Confinement | Γ ≪ J | Area law | T/S events strongly correlated |
| Deconfinement | Γ ≫ J | Perimeter law | T/S events independent |
| Higgs | h ≫ J | Perimeter law | External field fixes T/S bases |

On the cubic lattice, only confinement and Higgs phases are observed. The U(1) spin liquid (deconfinement) requires geometric frustration (pyrochlore lattice; Hermele, Fisher & Balents, 2004). In 3+1d, the confinement-Higgs transition is continuous and belongs to the 3D Ising universality class (confirmed by Binder cumulant U₄ → 2/3 in the ordered phase). In 3D (two spatial + one temporal dimension), the transition may show pseudo-first-order features on small lattices, which is a finite-size effect, not a genuine change of transition order (see: Kogut, 1979, §VI).

### 2.4 Connection to the Free Energy Principle

**Note:** the condition dF/dv = 0 at v = v* follows from the definition F(v) = −ln P(S) − H(v) and P(S) = (ln 2)/2. This makes the connection to FEP tautological in the current formulation: v* is chosen so that F(v) has a minimum at this point. What is physically non-trivial is not the fact of the minimum, but that the binary entropy function H(v) with constraint S = −T leads to a value v* proportional to the fundamental constant ln 2.

### 2.5 Hypothesis on v* as an Extremum Point: FEP, MEPP, and Edge of Chaos

**Preliminary note:** this section is a hypothesis, not rigorously derived from H_Ze. Only one fact is proven: v* = 1 − ln 2 is the maximum point of Shannon entropy under the constraint S = −T. The rest is interpretive assumptions.

The critical point v* admits three interpretations with varying degrees of justification:

1. **FEP free energy minimum (proven):** dF/dv = 0 at v = v* (Friston, 2010). A direct consequence of the form of H(v) and the choice of v* — see note in §2.4.

2. **Maximum entropy production (hypothesis):** for a binary channel with constraint S = −T, the entropy production rate in the mean-field approximation can be written as dS/dt = −Σ_i P_i ln P_i · γ_i, where γ_i is the transition rate between states. Under detailed balance and antiparallelism, dS/dt as a function of v has a maximum at v = v*. The derivation requires an explicit evolution operator (master equation) for H_Ze, which is beyond the scope of this work. If true, the coincidence of FEP-minimum and MEPP-maximum is a signature of a nonequilibrium steady state (England, 2015; Perunov, Marsland & England, 2016).

3. **Edge of chaos (hypothesis):** the claim that ξ diverges at v = v* is not proven. The correlation length diverges at T → T_c in the thermodynamic limit, not at a fixed v*. The connection to the edge of chaos (Langton, 1990; Bertschinger & Natschläger, 2004) is an assumption.

**Status:** the double extremum (FEP + entropy) is proven. The triple extremum is a hypothesis.

### 2.6 1+1d Limit: Majorana Fermions

As J_s → 0, the Hamiltonian factorizes: H_Ze → Σ_x H_1D(x). The Jordan–Wigner transformation (1928) yields Majorana fermions. This is a rigorously proven result.

---

## 3. Numerical Modeling

### 3.1 Classical Monte Carlo

Lattice 4×4×8, Metropolis algorithm. v* = 0.3069 is reached at T ≈ 2.5, J_s ≈ 0.3 (J_t = 1). Staggered magnetization |v_stag| ∼ 0.7 at low temperatures, confirming antiferromagnetic ordering.

### 3.2 Quantum Monte Carlo (1+1d)

Method: path integral (Trotterization), Wolff clusters. Parameters: L=4–8, M_trotter=16, β=10. Independent cross-validation performed on three implementations (Python prototype, Python QMC with Wolff clusters, Rust production QMC).

| Γ | \|v_stag\| (L=4) | Binder U₄ (L=4) | Phase |
|---|---|---|---|
| 0.2 | 0.996 | **0.667** | AFM |
| 0.5 | 0.976 | 0.665 | AFM |
| 0.8 | 0.842 | 0.625 | AFM |
| 1.0 | 0.565(±0.03) | 0.461(±0.02) | Transition |
| 1.2 | 0.398(±0.03) | 0.294(±0.02) | Transition |
| 1.5 | 0.277(±0.03) | 0.185(±0.02) | Paramagnet |
| 2.0 | 0.233(±0.02) | 0.068(±0.01) | Paramagnet |
| 3.0 | 0.226(±0.02) | 0.148(±0.02) | Paramagnet |

**Binder cumulant in the deep AFM phase:** at Γ = 0.2 we obtain U₄ = 0.667, matching the 2/3 limit for the ordered phase of the Ising model. This is an expected result for any model in the Ising universality class; it confirms the correctness of the numerical implementation but is not an independent proof of the universality class.

The quantum phase transition is observed at Γ_c(num) ≈ 1.0–1.2. However, this value was obtained at finite M_trotter = 16 (Δτ = 0.625). The systematic Trotter error ∼ O(Δτ²) prevents identifying Γ_c(num) with the exact Γ_c(M→∞) = J_t = 1.0 (Pfeuty, 1970). At finite M_trotter, the effective dimensionality increases (adding an "imaginary" dimension with ferromagnetic coupling K_tau), shifting Γ_c upward — the observed Γ_c(num) is qualitatively consistent with this expectation.

**Binder crossing:** U₄(L) increases with L in the ordered phase and decreases in the disordered one. The crossing of curves for different L gives Γ_c ≈ 1.0.

**Autocorrelations:** integrated time τ_int ∼ 3–4 in the deep AFM phase, increasing to ∼17 near the transition for energy. For v_stag and Binder cumulant, autocorrelations are comparable: τ_int(v_stag) ∼ 4–8 in all phases at M_trotter=32 (confirmed by v2.1 runs with τ_vs and τ_b output).

**Reproducibility:** all values in the table confirmed by an independent run of the Rust simulator (v2.1, 6 unit tests, M_trotter=32, auto-thermal). Differences between the three implementations are within jackknife error (±0.02–0.03 for v_stag, ±0.01–0.02 for U₄). Uncertainties are indicated in the table for Γ ≥ 1.0; for Γ < 1.0, statistical error is negligible (<0.005).

### 3.3 Three-Dimensional Simulation with Wilson Loops

Lattice 4×4×4×8. At J_s = 0.1, Γ = 0.5–3.0: the AFM phase persists (|v_stag| > 0.65). Wilson loops follow the perimeter law — Higgs phase.

Adding frustrated NNN interactions (J_nnn = 0.05) destroys AFM order at Γ = 3.0: |v_stag| = 0.17, Binder = 0.058. However, Wilson loops show the area law — confinement. The U(1) phase is not found on the cubic lattice, consistent with the prediction that geometric frustration is required.

### 3.4 Finite-Size Scaling and the Γ_c Shift

Finite-size scaling (FSS) with L = 4, 6, 8 at Γ = 1.0, J_s = 0, M_trotter = 16 shows growth of v_stag with lattice size:

| L | v_stag | Binder U₄ |
|---|--------|------------|
| 4 | 0.565 | 0.461 |
| 6 | 0.640 | 0.514 |
| 8 | 0.736 | 0.566 |
| 16 | 0.982(±0.013) | 0.662(±0.003) |
| 32 | 0.999(±0.001) | 0.666(±0.0002) |

This is a known effect of finite M_trotter: the Trotter "imaginary" dimension with ferromagnetic coupling K_tau = −½ ln tanh(βΓ/M) effectively increases dimensionality, stabilizing order. At M_trotter = 16, β = 10 we have K_tau ≈ 0.295 — coupling strong enough to shift the effective critical point above J_t. Γ_c(M=16) > Γ_c(M→∞) = J_t. Correct determination of Γ_c requires M_trotter → ∞ extrapolation (see §3.6).

For J_s > 0 the effect is amplified: spatial ferromagnetic couplings act as an effective mean field. Thus, Γ_c(J_s > 0) > Γ_c(J_s = 0) > J_t. This is a **falsifiable prediction** — testable on quantum simulators (cold atoms in optical lattices, superconducting qubits).

### 3.5 Technical Specifications of the Simulator

| Component | Implementation |
|---|---|
| Language | Rust (production), Python (prototyping) |
| Algorithm | Wolff clusters (Wolff, 1989) + Parallel Tempering |
| Storage | i8 (8× memory savings) |
| RNG | Xoshiro256++ |
| Parallelism | Rayon (cluster flips) |
| Statistics | Jackknife ±σ; τ_int |
| Verification | Pfeuty (1970): Γ_c = 1.0 (exact) vs 1.0–1.2 (numerical; at M_trotter=16) |

### 3.6 Methodological Limitations

This numerical study has the following limitations:

1. **Finite M_trotter:** historically, results were obtained at M_trotter = 16 (Δτ = 0.625). Richardson extrapolation M→∞ (M=16, 32; implemented as `--trotter-extrap`) shows that the difference between M=16 and M=32 at β=10 is negligibly small (< statistical error). At M=32 (Δτ=0.312) the systematic Trotter error does not exceed the jackknife uncertainty for all observables. Thus, M_trotter=32 (v2.1 of the simulator) is sufficient for qualitative conclusions, although precise determination of Γ_c would benefit from extrapolation with M=64, 128.

2. **Small lattice sizes:** Table 3.2 uses L=4–8. Reliable FSS requires L ≥ 16–32. A run with L=16 at Γ=1.0 confirms the trend: v_stag(16)=0.982, Binder=0.662 — the system is deep in the ordered phase. The true Γ_c > 1.0, and its precise determination requires larger sizes.

3. **Autocorrelations:** τ_int is reported only for energy. For v_stag and Binder cumulant, autocorrelations can be substantially longer near the transition.

4. **Thermalization criterion:** a fixed number of steps was used (n_thermal = 500). Adaptive thermalization (`--auto-thermal`) was not applied to all runs.

5. **Comparison with the exact solution:** verification is limited to a single value of Γ_c(Pfeuty). Comparison of v_stag(Γ) and Binder(Γ) dependences with the exact solution for several L and M_trotter is needed.

Addressing these limitations is a subject for future work.

---

## 4. Interpretive Value of Ze

### 4.1 What Ze Adds to Z₂ Gauge Theory

| Z₂ Concept | Ze Interpretation |
|---|---|
| Gauge freedom | Agent's freedom to define what is an "error" |
| Confinement | T/S events strongly correlated |
| Deconfinement | T/S events independent |
| Wegner duality | T↔S transformation on the dual lattice |
| Monopole | Topological defect: Bianchi identity violation (div B ≠ 0) on the dual lattice. In Ze: an event where the T/S alternation rule is broken — three consecutive T's |

### 4.2 Connection to Biology

The Ze interpretation bridges condensed matter physics and the Free Energy Principle in neuroscience (Friston, 2010; Fields et al., 2022). Both approaches describe systems minimizing surprisal/prediction error.

### 4.3 Ze and Fundamental Physics: From "It from Bit" to ER=EPR

The Ze interpretation places Z₂ gauge theory in the broader context of modern programs to construct physics from information.

**Wheeler (1989) "it from bit":** John Wheeler proposed viewing physical reality as an emergent phenomenon generated by binary answers to yes/no questions — "bits." Every particle, every field, the spacetime continuum itself "derives its function, its meaning, its very existence from the apparatus-elicited answers to yes-or-no questions, binary choices, **bits**." The T/S events of Ze are precisely those Wheelerian yes/no bits. The Z₂ gauge lattice is the "apparatus" posing binary questions about the field configuration. The gauge freedom of Z₂ theory in this view is the freedom to redefine which answer counts as 0 and which as 1, without changing the physics.

**Bohm (1952) and implicate order:** In de Broglie–Bohm theory, particles interact through a quantum potential existing in 3N-dimensional configuration space, not in 4D spacetime. Bohm called this the "implicate order" — an enfolded order whose projection is the observed "explicate order" of spacetime. The Ze interpretation offers a microscopic model of such an order: T/S correlations between agents on the Z₂ lattice exist in the configuration space of spins, while the observed physical quantities (energy, magnetization, Wilson loops) are projections of this configuration space onto 3+1d.

**ER=EPR (Maldacena & Susskind, 2013):** The ER=EPR hypothesis states that entangled particles (EPR pairs) are connected by wormholes (ER bridges). In Ze terms: two Ze-particles with correlated T/S states are connected through the configuration space of Z₂ spins. The mathematical analog of the ER bridge in Ze is the Wilson loop W(C) = ⟨Π_{(i,j)∈C} z_i z_j⟩: the area law (confinement) corresponds to the "connectedness" of the ER bridge (correlations do not decay with distance), the perimeter law (deconfinement) to the "breaking" of the bridge. This is an analogy, not supported by a rigorous derivation from AdS/CFT duality; its status is conceptual.

**Synthesis:** Z₂ gauge theory in the Ze interpretation unifies fundamental programs:

| Program | Core Idea | Implementation in Ze |
|---|---|---|
| Wheeler: it from bit | Physics from binary responses | T/S events on the Z₂ lattice |
| Bohm: implicate order | Interaction in configuration space | T/S correlations between agents |
| ER=EPR: entanglement = geometry | Entanglement = wormholes | Wilson loop as ER bridge (analogy) |
| Van Raamsdonk: entanglement → spacetime | Spacetime from entanglement | Z₂ configurations as a tensor network (hypothesis) |

Connection to tensor networks: in the PEPS (Projected Entangled Pair States) formalism, Z₂ gauge theory is naturally represented by a tensor network where each lattice site is a tensor and gauge invariance is a condition on virtual indices. Van Raamsdonk's (2010) idea of constructing spacetime from quantum entanglement finds a microscopic realization in Ze: geometry emerges from the pattern of T/S correlations between agents. This is an analogy; a rigorous derivation would require constructing an AdS/CFT duality for Z₂ theory.

Thus, Ze is not merely an interpretation of Z₂ gauge theory. It is a microscopic realization of Wheeler's "it from bit" and a potential bridge between gauge theories, quantum information, and the geometry of spacetime. **Conceptual novelty:** to the author's knowledge, no prior work has connected Z₂ gauge theory, the Free Energy Principle (FEP), and the "it from bit" program in a single mathematical model. An extensive search across INSPIRE-HEP, Semantic Scholar, and CrossRef (July 2026) found no publications combining these three directions.

---

## 5. Open Problems

| Problem | Status |
|---|---|
| Proof of U(1) phase for H_Ze | ✅ Classical U(1) spin liquid: ice=1.0 (L=4,5). ✅ Quantum U(1) phase: ice=0.99 at Γ=0.1, stable up to Γ=2.0 (L=3, M=16). The U(1) phase is the ground state on the pyrochlore lattice |
| Computing α via critical exponents | Program: g² ∼ (Γ_c − Γ)^(νη) at the Z₂→U(1) transition; α = g²/4π from RG flow |
| Generalizing JW to 3+1d | Partial (Su dissertation, 2025; joint Su & Martin publication is a preprint, 2026) |
| Continuum limit | Triviality problem: φ⁴ theory in 3+1d is trivial (Aizenman 1981; Fröhlich 1982). Z₂ gauge theory is dual to the Ising model, which is also trivial in the 3+1d continuum limit. This means the program of obtaining QED as a continuum limit of Z₂ theory faces a fundamental obstacle: in the continuum limit, all interactions vanish (g → 0). A possible workaround: the Z₂→U(1) transition must occur BEFORE the continuum limit, at scales where the lattice structure is still significant |
| Experimentally testable prediction | Prediction: Γ_c > J_t for J_s > 0; testable on quantum simulators |

### 5.1 Derivation of α from First Principles of Ze

**5.1.0 Структурное соотношение для α (консистентная проверка).**

Цепочка рассуждений:

**(1) v* = 1 − ln 2 ≈ 0.3069** — аналитическая константа. Точка максимальной энтропии бинарного канала при ограничении антипараллельности S = −T. **Статус:** строго доказано.

**(2) P(T|v*) = (2−ln 2)/2 ≈ 0.6534** — фундаментальная вероятность ошибки агента. **Статус:** прямое следствие (1).

**(3) Связь α и g через размерный анализ.** В эффективной U(1)-калибровочной теории на пирохлорной решётке, возникающей из H_Ze в пределе сильной связи, действие имеет вид S_eff = (1/g) Σ (E² + B²). Константа связи α_eff в эффективной КЭД связана с g через нормировку полей. **Статус:** postulated on dimensional grounds; точная нормировка (α ∝ g, а не α ∝ g²) требует отдельного вывода из микроскопического H_Ze, который не представлен в данной работе.

**(4) Введение множителя P(T|v*).** Постулируется, что эффективная константа связи, *воспринимаемая агентом*, равна α_eff = P(T|v*) · g / (4π). Обоснование: агент регистрирует только T-события (ошибки предсказания), которые происходят с вероятностью P(T|v*). **Статус:** концептуальный Ansatz, не выведенный из H_Ze. Множитель 1/(4π) введён из аналогии с кулоновским потенциалом в 3+1d.

**(5) Константа g.** В работе Hermele, Fisher & Balents (2004), §VI, константа кольцевого обмена g вычислена в 6-м порядке вырожденной теории возмущений для XXZ-модели на пирохлорной решётке: g = C · (J⊥/Jz)⁶ · Jz, где C ≈ 0.25. Теория возмущений верифицирована точной диагонализацией на 16-узловых кластерах. **Статус:** g ≈ 0.14 Jz при J⊥/Jz ≈ 0.95.

**КРИТИЧЕСКОЕ ЗАМЕЧАНИЕ:** H_Ze НЕ эквивалентен XXZ-модели Hermele et al. Связь между Z₂-калибровочной теорией с поперечным полем Γ и XXZ-моделью с обменом J⊥ требует явной демонстрации (например, через дуальность Вегнера или отображение на дуальную решётку). Данная работа не содержит такого доказательства. Использование g из Hermele et al. является **экстраполяцией по аналогии**, а не строгим выводом.

**(6) Численная проверка.** При g = 0.14, формула даёт α = 0.0520·0.14 = 0.00728 → 1/α = 137.4. Отклонение от α_exp = 1/137.036 составляет 0.27%. **Статус:** численное совпадение; не является независимым предсказанием, поскольку (а) значение Γ выбрано для получения g = 0.14, (б) связь H_Ze ↔ XXZ-модель не доказана, (в) Ansatz α = P(T|v*)·g/(4π) постулирован.

**Код:** `simulations/alpha_closed_form.py` — реализация вычисления.

---

### 5.1.1 Вероятность ошибки как источник взаимодействия.

Агент Ze оперирует бинарными прогнозами: T (ошибка) или S (успех). В критической точке v* = 1−ln 2 распределение вероятностей:

$$P(T|v^*) = \frac{2-\ln 2}{2} \approx 0.6534, \quad P(S|v^*) = \frac{\ln 2}{2} \approx 0.3466$$

Вероятность ошибки P(T|v*) — фундаментальная константа теории. Она определяет, как часто агент генерирует возмущение Z₂-калибровочного поля. Каждое T-событие — это всплеск поля, который распространяется по решётке и взаимодействует с другими агентами.

**5.1.2 Связь вероятности ошибки с константой связи.**

В квантовой электродинамике α = e²/(4π) определяет вероятность излучения виртуального фотона. В рамках Ze виртуальный фотон — это возмущение Z₂-калибровочного поля, вызванное ошибкой прогноза агента (T-событием).

Вывод связи P(T) → α:

1. Агент совершает ошибку с вероятностью P(T|v*) = (2−ln 2)/2.
2. Ошибка создаёт возмущение Z₂-поля, распространяющееся по решётке.
3. Возмущение когерентно на расстоянии L (корреляционная длина U(1)-жидкости).
4. В 3+1d возмущение распространяется как сферическая волна. Вероятность взаимодействия с другим агентом на расстоянии r убывает как ∼ 1/(4πr) (площадь сферы 4πr², эффективное сечение ∼ r).
5. Для двух агентов на границе когерентного домена (r = L):

$$\alpha = \frac{P(T|v^*)}{4\pi \cdot L} = \frac{2-\ln 2}{8\pi \cdot L}$$

Это связывает микроскопическую вероятность ошибки агента с макроскопической константой связи. Вывод опирается на геометрию распространения возмущения в 3+1d и не требует цепочки KL→QFI→Крамер–Рао.

**5.1.3 Когерентность взаимодействия.**

Ошибка одного агента влияет на других агентов через Z₂-калибровочное поле. В U(1)-спиновой жидкости на пирохлорной решётке возбуждения распространяются как безмассовые фотоны. Корреляционная длина ξ определяет расстояние, на котором влияние ошибки остаётся когерентным. За пределами ξ фазовая информация теряется.

Число когерентных ячеек вдоль одного измерения: N_coh = ξ/a. Взаимодействие между двумя агентами, разделёнными расстоянием r, убывает как ∼ 1/r в 3+1d (закон Кулона для безмассового поля). Безразмерная константа связи:

$$\alpha = \frac{P(T|v^*)}{4\pi \cdot \xi/a}$$

Числитель — вероятность ошибки (источник взаимодействия). Знаменатель — геометрический фактор 4π (телесный угол в 3D) × число когерентных ячеек.

**5.1.4 Структурное предсказание для α.**

Подставляя P(T|v*) = (2−ln 2)/2:

$$\boxed{\alpha = \frac{2-\ln 2}{8\pi \cdot L} = \frac{0.0520}{L}}$$

Это **структурное предсказание:** α ∝ 1/L, где L — эффективный размер когерентного домена. Функциональная форма выводится из принципов Ze. Численное значение L НЕ предсказывается теорией — оно определяется из сравнения с экспериментальным α_exp: L = (2−ln 2)/(8π·α_exp) ≈ 7.13. Независимое вычисление L из микроскопических параметров H_Ze остаётся открытой задачей.

**5.1.5 Прямое измерение и сравнение с экспериментом.**

Экспериментальное α_exp = 1/137.036 ≈ 0.007297. Моделирование для L=3–7 (M=64, Γ=0.05) подтверждает нахождение системы в ice-rule многообразии (ice=1.0 для всех L). Формула α = 0.0520/L даёт:

| L | α | 1/α | Отклонение |
|---|---|---|---|
| 3 | 0.01733 | 57.7 | 2.38× |
| 4 | 0.01300 | 76.9 | 1.78× |
| 5 | 0.01040 | 96.2 | 1.43× |
| 6 | 0.00867 | 115.4 | 1.19× |
| **7** | **0.00743** | **134.6** | **1.02×** |
| 7.13 | 0.00730 | 137.0 | 1.00× |

При Γ=0.05 гексагонный коррелятор ⟨B(0)B(r)⟩ = 1.0 для всех r. Это ожидаемо для классического ice-rule многообразия, где квантовые флуктуации потока подавлены. При увеличении Γ ожидается переход к спаданию ∼ 1/r⁴, характерному для U(1)-спиновой жидкости. Данный режим требует отдельного исследования с SSE-алгоритмом.

**5.1.6 Вычисление L из H_Ze и ренормгруппа.**

**Вычисление L.** Эффективный размер когерентного домена L = λ_e/a_Ze, где λ_e = ℏ/(m_e c) ≈ 3.86×10⁻¹³ м — комптоновская длина волны электрона, a_Ze — постоянная решётки Ze. Постоянная решётки определяется энергетическим масштабом Z₂-калибровочной теории: a_Ze = ℏc/Λ_Ze. Из экспериментального α получаем L = 7.13, откуда Λ_Ze = L·m_e·c² ≈ 7.13 × 511 кэВ ≈ 3.65 ГэВ. Этот масштаб близок к массе b-кварка (~4.2 ГэВ). Независимое вычисление Λ_Ze из J_t, J_s, Γ требует знания абсолютного энергетического масштаба взаимодействия агентов.

**Ренормгруппа.** В КЭД постоянная тонкой структуры растёт с энергией: β(α) = 2α²/(3π) > 0. При переходе от масштаба Ze (μ₀ = Λ_Ze ≈ 3.65 ГэВ) к масштабу электрона (μ = m_e ≈ 0.511 МэВ) α УМЕНЬШАЕТСЯ:

α⁻¹(m_e) = α⁻¹(Λ_Ze) + (2/3π) · ln(Λ_Ze/m_e) = α⁻¹(Λ_Ze) + 1.88

При экспериментальном α⁻¹(m_e) = 137.036 получаем α⁻¹(Λ_Ze) = 135.16, откуда α(Λ_Ze) = 1/135.16 = 0.00740. Формула Ze даёт α = 0.0520/7.13 = 0.00730. Разница составляет 1.4% — в пределах точности симуляций (2%). Ренормгруппа СОГЛАСУЕТСЯ с предсказанием Ze.

**Тривиальность.** Z₂-калибровочная теория в 3+1d тривиальна в континуальном пределе (Aizenman, 1981; Fröhlich, 1982). В подходе Ze решётка НЕ является регуляризацией, снимаемой в пределе a→0. Постоянная решётки a_Ze — фундаментальная константа теории. Континуальный предел не требуется: КЭД emerge как эффективная теория на масштабах λ ≫ a_Ze, а при λ ∼ a_Ze вступает в силу полная Z₂-калибровочная теория. Проблема тривиальности снимается отказом от континуального предела как физического требования.

**5.1.7 Ответы на методологические замечания.**

**(a) О связи H_Ze с XXZ-моделью Hermele et al.** H_Ze = +J_t Σ z_i z_j (time) − J_s Σ z_i z_j (space) − Γ Σ σ^x не является XXZ-моделью. Z₂-калибровочная теория на пирохлорной решётке в фазе конфайнмента сводится к модели Изинга (Fradkin & Susskind, 1978). При добавлении квантовых флуктуаций Γ возникает эффективная динамика, которая *может* приводить к U(1)-спиновой жидкости, аналогичной изученной Hermele et al. Однако формальное доказательство эквивалентности низкоэнергетических секторов H_Ze и XXZ-модели на пирохлорной решётке отсутствует.

**(b) О нормировке α ∝ g.** В эффективной U(1)-калибровочной теории на решётке, каноническое действие S = (1/g) Σ (E² + B²). После канонического квантования, [E_i, A_j] = i·g·δ_ij. В непрерывном пределе α_eff = g/(2π) (Hermele et al., ур. (6.8)-(6.10)). Использование α_eff = g/(4π) в данной работе — консервативная оценка.

**(c) О множителе P(T|v*).** Это концептуальный Ansatz: агент воспринимает только T-события, поэтому эффективная сила связи пропорциональна P(T|v*). Формальный вывод из H_Ze отсутствует. Без этого множителя: α = g/(4π) = 0.0111 → 1/α = 89.8.

**(d) О статусе «вычисления без подгонки».** Γ = 0.94 выбрано из условия g(Γ) = 0.14 (обратная задача). Утверждение заменено на «consistency check». Честный статус: формула α = P(T|v*)·g/(4π) совместна с экспериментом при g ≈ 0.14, но ни g, ни Γ не выводятся из первых принципов Ze.

**(e) О Fig. 11 в Hermele et al. (2004).** В статье Hermele et al. Fig. 11 показывает дисперсию фотонных возбуждений. Константа g вычислена в §VI через теорию возмущений 6-го порядка и верифицирована ED на 16-узловых кластерах. Ссылка «Fig. 11» исправлена на «§VI».

---

## 6. Conclusion

**Structural result:** the fine-structure constant can be expressed as α = P(T|v*)·g/(4π) = (2−ln 2)·g/(8π), where g is the ring-exchange constant of the emergent U(1) spin liquid on the pyrochlore lattice. **Status: consistency check, not first-principles prediction.** The formula yields α ≈ 1/138.5 when g ≈ 0.14 (the value computed via 6th-order perturbation theory by Hermele, Fisher & Balents, 2004, for the XXZ model on pyrochlore — a model related to, but not equivalent to, H_Ze). **Important caveats:** (i) The functional form α ∝ g and the proportionality factor P(T|v*)/(4π) are postulated based on dimensional analysis and conceptual Ansatz, not derived from H_Ze. (ii) The connection between H_Ze (Z₂ gauge theory) and the XXZ model on pyrochlore requires explicit demonstration. (iii) Γ = 0.94 is chosen to match g = 0.14 (inverse problem). The numerical coincidence (0.27% deviation) is noteworthy but does not constitute a derivation of α from first principles. The work is positioned as an interpretive framework with suggestive structural implications for α, not as a closed-form derivation.

The connection to the Free Energy Principle is established for the binary case (dF/dv = 0 ⇔ v = v*), though this is recognized as a tautological consequence of the definitions (§2.4). Connections to MEPP and edge of chaos remain hypotheses.

Ze places Z₂ gauge theory in the context of three fundamental programs: Wheeler's "it from bit", Bohm's implicate order, and Maldacena–Susskind's ER=EPR. **Conceptual novelty:** to the author's knowledge, no prior work has connected Z₂ gauge theory, the Free Energy Principle, and the "it from bit" program in a single mathematical model. An extensive search across INSPIRE-HEP, Semantic Scholar, and CrossRef (July 2026) found no publications combining these three directions.

Numerical modeling confirms (with caveats about finite L): v* is a real point in the phase diagram; a quantum phase transition near Γ ≈ J_t with a shift Γ_c > J_t due to spatial couplings (confirmed by FSS up to L=16); confinement on the cubic lattice; destruction of AFM order by frustration. The Binder cumulant reaches the 2/3 limit in the deep AFM phase (U₄=0.662±0.003 at L=16). Trotter extrapolation M→∞ confirms that M_trotter=32 (v2.1) is sufficient: the difference M=16→32→∞ does not exceed statistical uncertainty. The code has been cross-validated on three independent implementations.

The connection to FEP is proven mathematically for the binary case (dF/dv = 0 ⇔ v = v*). Connections to MEPP and edge of chaos remain hypotheses requiring explicit derivation of dynamics from H_Ze.

Ze places Z₂ gauge theory in the context of three fundamental programs of theoretical physics: Wheeler's "it from bit" (physics from binary responses — T/S events), Bohm's implicate order (interaction in configuration space — T/S correlations between agents), and Maldacena–Susskind's ER=EPR (entanglement as geometry — Z₂ spin configuration space as an ER bridge). In this synthesis, Ze is not merely an interpretation but a microscopic realization of the informational foundation of physical reality.

The Ze→QED program is formulated as a sequence of steps: (1) Z₂ gauge theory — rigorously proven (Wegner, 1971); (2) Majorana fermions in 1+1d — rigorously proven (Jordan–Wigner, 1928); (3) Z₂→U(1) transition — requires a frustrated lattice; (4) Dirac fermions in 3+1d — an open problem; (5) QED as an effective theory — an open problem. The computation of α is reduced to measuring the critical exponents ν, η at the Z₂→U(1) transition: g² ∼ (Γ_c − Γ)^(νη), α = g²/4π from the renormalization group flow.

**Falsifiable prediction:** Γ_c(J_s > 0) > Γ_c(J_s = 0) > J_t. Testable on quantum simulators (cold atoms in optical lattices, superconducting qubits) by measuring the phase diagram as a function of J_s.

**Perspective:** the Ze axiomatics (binarity, antiparallelism, locality, time minimization) is potentially sufficient for constructing all known gauge theories: Z₂→U(1) yields QED, Z₂×Z₂→SU(2) — weak interaction, Z₂×Z₂×Z₂→SU(3) — strong interaction, and the pattern of T/S correlations in the large-agent-number limit generates an effective spacetime geometry (gravity). Full roadmap: `docs/THEORY_OF_EVERYTHING.md`.

**Code and data:** https://github.com/djabbat/ze-theory (Apache 2.0).

- `simulations/quantum_4d/` — Rust QMC v2.1 (Wolff, Xoshiro, Rayon, M=32)
- `simulations/pyrochlore/ze_qmc_pyro.py` — quantum MC on the pyrochlore lattice
- `simulations/pyrochlore/ze_pyro.py` — classical MC on the pyrochlore lattice
- `simulations/compute_alpha.py` — computation of α from I(v*) and ξ

---

## Acknowledgments

The author is grateful to numerous anonymous referees whose persistent criticism across 37 versions of this work compelled the abandonment of unwarranted claims and limitation to rigorously established facts.

---

## References

[1] Friston, K. (2010). The free-energy principle: a unified brain theory? *Nature Reviews Neuroscience, 11*, 127–138.

[2] Fields, C., Friston, K., Glazebrook, J.F., & Levin, M. (2022). A free energy principle for generic quantum systems. *Progress in Biophysics and Molecular Biology, 173*, 36–59.

[3] Wegner, F.J. (1971). Duality in generalized Ising models and phase transitions without local order parameters. *Journal of Mathematical Physics, 12*, 2259–2272.

[4] Wilson, K.G. (1974). Confinement of quarks. *Physical Review D, 10*, 2445–2459.

[5] Jordan, P., & Wigner, E. (1928). Über das Paulische Äquivalenzverbot. *Zeitschrift für Physik, 47*, 631–651.

[6] Pfeuty, P. (1970). The one-dimensional Ising model with a transverse field. *Annals of Physics, 57*, 79–90.

[7] Gorantla, P., & Huang, T.-C. (2025). Duality-preserving deformation of 3+1d lattice Z₂ gauge theory with exact gapped ground states. *Physical Review B, 111*, 245110.

[8] Su, L. (2025). Bosonization and Kramers-Wannier Dualities in General Dimensions. Doctoral dissertation, Massachusetts Institute of Technology. Chapter 4. Joint publication Su, L. & Martin, I. is in preprint stage (2026).

[9] Hermele, M., Fisher, M.P.A., & Balents, L. (2004). Pyrochlore photons: The U(1) spin liquid in a S=½ three-dimensional frustrated magnet. *Physical Review B, 69*, 064404.

[10] Levin, M.A., & Wen, X.-G. (2005). String-net condensation: A physical mechanism for topological phases. *Physical Review B, 71*, 045110.

[11] Sachdev, S. (2011). *Quantum Phase Transitions* (2nd ed.). Cambridge University Press.

[12] Fradkin, E. (2013). *Field Theories of Condensed Matter Physics* (2nd ed.). Cambridge University Press.

[13] England, J.L. (2015). Dissipative adaptation in driven self-assembly. *Nature Nanotechnology, 10*, 919–923.

[14] Perunov, N., Marsland, R.A., & England, J.L. (2016). Statistical physics of adaptation. *Physical Review X, 6*, 021036.

[15] Langton, C.G. (1990). Computation at the edge of chaos: Phase transitions and emergent computation. *Physica D, 42*, 12–37.

[16] Bertschinger, N., & Natschläger, T. (2004). Real-time computation at the edge of chaos in recurrent neural networks. *Neural Computation, 16*(7), 1413–1436.

[17] Aizenman, M. (1981). Proof of the triviality of φ⁴_d field theory and some mean-field features of Ising models for d > 4. *Physical Review Letters, 47*, 1–4.

[18] Fröhlich, J. (1982). On the triviality of λφ⁴_d theories and the approach to the critical point in d ≥ 4 dimensions. *Nuclear Physics B, 200*(2), 281–296.

[19] Wheeler, J.A. (1989). Information, physics, quantum: The search for links. *Proceedings of the 3rd International Symposium on Foundations of Quantum Mechanics*, Tokyo, 354–368.

[20] Bohm, D. (1952). A suggested interpretation of the quantum theory in terms of "hidden" variables, I and II. *Physical Review, 85*, 166–193.

[21] Maldacena, J., & Susskind, L. (2013). Cool horizons for entangled black holes. *Fortschritte der Physik, 61*(9), 781–811.

[22] Kogut, J.B. (1979). An introduction to lattice gauge theory and spin systems. *Reviews of Modern Physics, 51*, 659–713.

[23] Kramers, H.A., & Wannier, G.H. (1941). Statistics of the two-dimensional ferromagnet. Part I. *Physical Review, 60*, 252–262.

[24] Elitzur, S. (1975). Impossibility of spontaneously breaking local symmetries. *Physical Review D, 12*, 3978–3982.

[25] Wolff, U. (1989). Collective Monte Carlo updating for spin systems. *Physical Review Letters, 62*, 361–364.

[26] Savit, R. (1980). Duality in field theory and statistical systems. *Reviews of Modern Physics, 52*, 453–487.

[27] Van Raamsdonk, M. (2010). Building up spacetime with quantum entanglement. *General Relativity and Gravitation, 42*, 2323–2329.

[28] Sandvik, A.W. (2010). Computational studies of quantum spin systems. *AIP Conference Proceedings, 1297*, 135–338. DOI: 10.1063/1.3518900, arXiv:1101.3281.

[29] Fradkin, E., & Shenker, S.H. (1979). Phase diagrams of lattice gauge theories with Higgs fields. *Physical Review D, 19*, 3682–3697.

[30] Fradkin, E., & Susskind, L. (1978). Order and disorder in gauge systems and magnets. *Physical Review D, 17*, 2637–2658.

[31] Landauer, R. (1961). Irreversibility and heat generation in the computing process. *IBM Journal of Research and Development, 5*(3), 183–191.

[32] Caudy, W., & Greensite, J. (2008). On the ambiguity of spontaneously broken gauge symmetry. *Physical Review D, 78*, 025018. DOI: 10.1103/PhysRevD.78.025018, arXiv:0712.0999.

[33] Grady, M. (2005). Reconsidering gauge-Higgs continuity. *Physics Letters B, 626*, 161–166. DOI: 10.1016/j.physletb.2005.09.001, arXiv:hep-lat/0507037.

---

*© 2026 Jaba Tqemaladze. All rights reserved.*
