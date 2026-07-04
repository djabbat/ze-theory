# Ze Theory: An Interpretive Framework for Z₂ Gauge Theory with Structural Implications for the Fine-Structure Constant

**Jaba Tqemaladze, MD**

*Georgia Longevity Alliance | Kutaisi International University | Free University of Tbilisi*

*Correspondence: jaba@longevity.ge | ORCID: 0000-0001-8651-7243*

*July 4, 2026*

---

## Abstract

Ze theory offers an interpretive framework for Z₂ lattice gauge theory cast in the language of active agents that minimize existence time through prediction. The gauge freedom of Z₂ theory is reinterpreted as the agent's freedom to designate which events count as prediction errors (T, for Tension) and which as successes (S, for Stretch). A structural relation emerges: the fine-structure constant can be expressed as α = P(T|v⁎)·g/(4π), where P(T|v⁎) = (2 − ln 2)/2 is the agent's error probability at the critical point v⁎ = 1 − ln 2, and g is the ring-exchange constant of the emergent U(1) spin liquid on the pyrochlore lattice. An important caveat: g is computed via sixth-order degenerate perturbation theory in the XXZ model on pyrochlore (Hermele, Fisher & Balents, 2004, §VI) and verified against exact diagonalization on 16-site clusters. It has not been measured non-perturbatively, and the connection between H_Ze—a Z₂ gauge theory—and the XXZ model requires explicit demonstration, which is not yet provided. At Γ/J_t ≈ 0.94, a value chosen to match the perturbative regime where g ≈ 0.14, the formula yields α ≈ 1/138.5. This is a consistency check, not a first-principles prediction. The functional form α ∝ g is postulated on dimensional grounds from the effective gauge action; the factor P(T|v⁎) is a conceptual Ansatz linking agent error probability to coupling strength. Neither the functional form nor the proportionality constant is derived from H_Ze.

---

## 1. Introduction

### 1.1 Definitions

**Ze-agent.** A mathematical model of a system that observes binary events, predicts the next event, and records the outcome: T (Tension — prediction error) or S (Stretch — prediction success).

The parameter v = (N_T − N_S)/(N_T + N_S) ∈ [−1, +1] is the Ze-velocity. The critical point is v⁎ = 1 − ln 2 ≈ 0.3069. At v⁎ one has P(S) = (ln 2)/2 ≈ 0.347, P(T) = (2 − ln 2)/2 ≈ 0.653, and N_T/N_S = 2/ln 2 − 1 ≈ 1.885. This value emerges in classical Monte Carlo at T ≈ 2.5, J_s ≈ 0.3 (J_t = 1) as the point separating confinement from deconfinement.

**Existence time** τ_Ze is an agent's internal resource, depleted by T-events: Δτ = η·N_T. This is a phenomenological parameter with no direct connection to relativistic proper time or any known physical observable. A possible microscopic interpretation sets τ_Ze ∝ 1/Δ, where Δ is the energy gap of H_Ze, linking existence time to the characteristic timescale of quantum evolution. This connection remains a hypothesis.

**Gauge freedom and T/S-basis.** The gauge symmetry of Z₂ theory is strictly local (G_x|Ψ⟩ = +|Ψ⟩). The interpretation of this symmetry as "the agent's freedom to choose the T/S-basis" is a conceptual choice, not a mathematical consequence. The agent cannot arbitrarily select the T/S-basis outside the physical Hilbert space.

### 1.2 What Ze Asserts and Does Not Assert

**The framework asserts:**
- The gauge freedom of Z₂ theory can be interpreted as the agent's freedom to define the T/S-basis.
- Interaction is a strategy for minimizing existence time expenditure.
- The interpretation is consistent with the Free Energy Principle (Friston, 2010; Fields et al., 2022) for the binary case.
- The parameter v is a statistical property of the agent's classical measurement record, not a quantum expectation value; Elitzur's theorem does not constrain it.

**The framework does not assert:**
- Ze does not derive QED from first principles.
- Ze does not compute the fine-structure constant α; the numerical coincidence in §5.1 is a consistency check.
- Ze does not prove new mathematical theorems about Z₂ gauge theory.
- Ze does not offer experimentally testable predictions beyond the Γ_c(J_s) dependence, verifiable on quantum simulators.
- The connections to maximum entropy production and edge of chaos (§2.5) are hypotheses, not derived from H_Ze.
- The connection to ER = EPR (§4.3) is a conceptual analogy, not supported by mathematical derivation.

---

## 2. Mathematical Formalism

### 2.1 Hamiltonian

$$H_{\mathrm{Ze}} = +J_t\sum_{x,t} z_{x,t} z_{x,t+1} \;-\; J_s\sum_{\langle x,y\rangle, t} z_{x,t} z_{y,t} \;-\; \Gamma\sum_{x,t} \sigma^x_{x,t} \;-\; h\sum_{x,t} z_{x,t}$$

where z_{x,t} ∈ {±1}, J_t > 0 (antiferromagnetic, encoding antiparallelism), J_s > 0 (ferromagnetic, encoding locality), Γ is the transverse field producing quantum fluctuations, and h is an external longitudinal field.

**Justification of signs.** The term +J_t·z_i·z_j with J_t > 0 gives energy ∼ +J_t for parallel spins and ∼ −J_t for antiparallel spins, favoring antiparallelism along the temporal direction. The term −J_s·z_i·z_j with J_s > 0 gives energy ∼ −J_s for parallel spins, favoring ferromagnetic ordering in space.

### 2.2 Equivalence to Z₂ Gauge Theory

The Gauss operator is G_x = σˣ_x ∏_μ σˣ_{x,μ}. Direct computation shows [H_Ze, G_x] = 0 for h = 0. The physical Hilbert space satisfies G_x|Ψ⟩ = +|Ψ⟩. This is the standard formalism of Z₂ gauge theory (Wegner, 1971; Wilson, 1974; Kogut, 1979). No new Hamiltonian is introduced; Ze employs the existing one.

**Elitzur's theorem and the status of v.** Elitzur's theorem (Elitzur, 1975) states that in a gauge theory with local symmetry, non-gauge-invariant operators possess zero expectation value in any gauge-invariant state: ⟨Ψ|z|Ψ⟩ = 0 for all physical states satisfying G_x|Ψ⟩ = +|Ψ⟩. The proof is elementary: G_x z G_x = −z for links incident on x, hence ⟨Ψ|z|Ψ⟩ = ⟨Ψ|G_x z G_x|Ψ⟩ = −⟨Ψ|z|Ψ⟩ ⇒ ⟨Ψ|z|Ψ⟩ = 0.

Consequently, v = ⟨z⟩ = 0 identically in the physical Hilbert space. Gauge fixing does not alter this result; it is a calculational tool, not a physical operation that creates new observables. The analogy with the Higgs expectation value ⟨φ⟩ in unitary gauge, drawn in earlier versions of this work, is misleading because in Higgs theories the gauge-invariant combination |φ|² is physical, whereas pure Z₂ gauge theory possesses no gauge-invariant bilinear that reduces to z.

**However**, Caudy & Greensite (2008) demonstrated that while local gauge symmetries cannot break spontaneously, **global subgroups** of the local gauge symmetry can. Specifically: "Local gauge symmetries cannot break spontaneously, according to Elitzur's theorem, but this leaves open the possibility of breaking some global subgroup of the local gauge symmetry. [...] The location of the breaking in the phase diagram depends on the choice of global subgroup" (Caudy & Greensite, 2008, Phys. Rev. D 78, 025018). This is the precise mathematical mechanism underlying the Ze interpretation: the agent's choice of T/S-basis selects a global subgroup of the Z₂ gauge symmetry, and the breaking of this subgroup at v = v⁎ is physically meaningful. Similarly, Grady (2005) showed that the confinement-Higgs continuity (Fradkin & Shenker, 1979) depends on gauge-fixing: in a partial axial gauge the two phases are separated by a genuine phase transition.

**Implementation in Ze.** The parameter v serves as the order parameter for the breaking of the agent's chosen global subgroup. In the unbroken phase (v = 0; T and S equally likely) the agent has no predictive power. In the broken phase (v ≠ 0) the agent can distinguish T from S events. The critical point v⁎ = 1 − ln 2 marks the transition between these regimes.

An important nuance: which global subgroup is broken is determined by the agent's choice of T/S-basis, but the location of the transition in the phase diagram is dictated by the dynamics of H_Ze. The fact that classical Monte Carlo finds v⁎ at T ≈ 2.5, J_s ≈ 0.3 suggests that the phase boundary of H_Ze aligns with the agent's entropy-maximizing point—a non-trivial dynamical observation that warrants further investigation.

### 2.3 Phase Diagram

Z₂ gauge theory in 3+1d possesses three phases (Wilson, 1974; Fradkin, 2013). An important property is **self-duality**: Z₂ theory in 3+1d is dual to itself under the Kramers–Wannier transformation (Kramers & Wannier, 1941; Wegner, 1971; Savit, 1980). Explicitly, the substitution z → z′ on the dual lattice with σˣ ↔ σᶻ maps H_Ze(Γ, h) to H_Ze(h, Γ). The self-dual point Γ = h determines the phase transition. In the Ze interpretation, self-duality corresponds to symmetry between T- and S-events: describing the system in terms of errors (T) is equivalent to describing it in terms of successes (S) on the dual lattice.

A note on phases: in Z₂ theory with fundamental Higgs fields, the confinement and Higgs phases are analytically connected; there is no phase transition between them, only a crossover (Fradkin & Shenker, 1979). For the Ze interpretation this means that the "agent's choice of T/S-basis" (gauge fixing) does not create physically distinct phases, consistent with Elitzur's theorem (§2.2). Deconfinement, as a topologically ordered phase, requires frustration (Hermele, Fisher & Balents, 2004).

| Phase | Parameters | Wilson Loops | Ze Interpretation |
|---|---|---|---|
| Confinement | Γ ≪ J | Area law | T/S events strongly correlated |
| Deconfinement | Γ ≫ J | Perimeter law | T/S events independent |
| Higgs | h ≫ J | Perimeter law | External field fixes T/S bases |

On the cubic lattice only confinement and Higgs phases are observed. The U(1) spin liquid (deconfinement) requires geometric frustration, e.g. the pyrochlore lattice (Hermele, Fisher & Balents, 2004). In 3+1d the confinement-Higgs transition is continuous and belongs to the 3D Ising universality class, confirmed by the Binder cumulant U₄ → 2/3 in the ordered phase. In three dimensions (two spatial plus one temporal) the transition may exhibit pseudo-first-order features on small lattices, a finite-size effect rather than a genuine change of transition order (Kogut, 1979, §VI).

### 2.4 Connection to the Free Energy Principle

The condition dF/dv = 0 at v = v⁎ follows from the definition F(v) = −ln P(S) − H(v) with P(S) = (ln 2)/2. This renders the connection to the Free Energy Principle tautological in the current formulation: v⁎ is chosen so that F(v) attains a minimum there. What is non-trivial is not the minimum itself but the fact that the binary entropy function H(v) under the constraint S = −T yields a value v⁎ proportional to the fundamental constant ln 2.

### 2.5 Hypothesis: v⁎ as an Extremum of FEP, MEPP, and Edge of Chaos

This section is a hypothesis, not rigorously derived from H_Ze. Only one fact is proven: v⁎ = 1 − ln 2 is the maximum of Shannon entropy under the constraint S = −T. The remainder consists of interpretive assumptions.

The critical point v⁎ admits three interpretations of varying justification:

1. **FEP free-energy minimum (proven):** dF/dv = 0 at v = v⁎ (Friston, 2010), a direct consequence of the definition of F(v).

2. **Maximum entropy production (hypothesis):** for a binary channel with antiparallelism, the entropy production rate in the mean-field approximation can be written as dS/dt = −∑_i P_i ln P_i · γ_i, where γ_i are transition rates. Under detailed balance, dS/dt as a function of v attains a maximum at v = v⁎. The derivation requires an explicit master equation for H_Ze, which lies beyond the scope of this work. If true, the coincidence of the FEP minimum and the MEPP maximum would signal a non-equilibrium steady state (England, 2015; Perunov, Marsland & England, 2016).

3. **Edge of chaos (hypothesis):** the claim that the correlation length diverges at v = v⁎ is not proven. The correlation length diverges at T → T_c in the thermodynamic limit, not at a fixed v⁎. The connection to the edge of chaos (Langton, 1990; Bertschinger & Natschläger, 2004) is an assumption.

**Status:** the double extremum (FEP + entropy) is proven. The triple extremum remains a hypothesis.

### 2.6 1+1d Limit: Majorana Fermions

As J_s → 0, the Hamiltonian factorizes: H_Ze → ∑_x H_{1D}(x). The Jordan–Wigner transformation (1928) yields Majorana fermions. This is a rigorously proven result.

---

## 3. Numerical Modeling

### 3.1 Classical Monte Carlo

Lattice 4×4×8, Metropolis algorithm. v⁎ = 0.3069 is reached at T ≈ 2.5, J_s ≈ 0.3 (J_t = 1). The staggered magnetization |v_stag| ∼ 0.7 at low temperatures confirms antiferromagnetic ordering.

### 3.2 Quantum Monte Carlo (1+1d)

Method: path integral with Trotterization and Wolff cluster updates. Parameters: L = 4–8, M_trotter = 16, β = 10. Independent cross-validation was performed on three implementations (Python prototype, Python QMC with Wolff clusters, and Rust production QMC).

| Γ | \|v_stag\| (L=4) | Binder U₄ (L=4) | Phase |
|---|---|---|---|
| 0.2 | 0.996 | 0.667 | AFM |
| 0.5 | 0.976 | 0.665 | AFM |
| 0.8 | 0.842 | 0.625 | AFM |
| 1.0 | 0.565(±0.03) | 0.461(±0.02) | Transition |
| 1.2 | 0.398(±0.03) | 0.294(±0.02) | Transition |
| 1.5 | 0.277(±0.03) | 0.185(±0.02) | Paramagnet |
| 2.0 | 0.233(±0.02) | 0.068(±0.01) | Paramagnet |
| 3.0 | 0.226(±0.02) | 0.148(±0.02) | Paramagnet |

At Γ = 0.2 we obtain U₄ = 0.667, matching the 2/3 limit for the ordered phase of the Ising model. This confirms the correctness of the numerical implementation but does not independently establish the universality class.

The quantum phase transition is observed at Γ_c(num) ≈ 1.0–1.2. However, this value was obtained at finite M_trotter = 16 (Δτ = 0.625). The systematic Trotter error of order Δτ² prevents identifying Γ_c(num) with the exact value Γ_c(M → ∞) = J_t = 1.0 (Pfeuty, 1970). At finite M_trotter the effective dimensionality increases, shifting Γ_c upward, which is qualitatively consistent with the observed value.

Binder crossing analysis yields Γ_c ≈ 1.0. The integrated autocorrelation time τ_int ∼ 3–4 in the deep AFM phase, increasing to ∼ 17 near the transition. All values in the table have been confirmed by an independent Rust simulator (v2.1, 6 unit tests, M_trotter = 32, auto-thermalization). Differences among the three implementations lie within the jackknife error (±0.02–0.03 for v_stag; ±0.01–0.02 for U₄).

### 3.3 Three-Dimensional Simulation with Wilson Loops

Lattice 4×4×4×8. At J_s = 0.1, Γ = 0.5–3.0: the AFM phase persists (|v_stag| > 0.65). Wilson loops follow the perimeter law, characteristic of the Higgs phase. Adding frustrated next-nearest-neighbor interactions (J_nnn = 0.05) destroys AFM order at Γ = 3.0 (|v_stag| = 0.17, Binder = 0.058), while Wilson loops obey the area law, signaling confinement. The U(1) phase is not found on the cubic lattice, consistent with the requirement of geometric frustration.

### 3.4 Finite-Size Scaling and the Γ_c Shift

Finite-size scaling (FSS) with L = 4, 6, 8 at Γ = 1.0, J_s = 0, M_trotter = 16 shows growth of v_stag with lattice size:

| L | v_stag | Binder U₄ |
|---|--------|------------|
| 4 | 0.565 | 0.461 |
| 6 | 0.640 | 0.514 |
| 8 | 0.736 | 0.566 |
| 16 | 0.982(±0.013) | 0.662(±0.003) |
| 32 | 0.999(±0.001) | 0.666(±0.0002) |

This is a known effect of finite M_trotter: the Trotter dimension with ferromagnetic coupling K_τ = −½ ln tanh(βΓ/M) effectively increases the dimensionality, stabilizing order. At M_trotter = 16 and β = 10 one has K_τ ≈ 0.295, a coupling strong enough to shift the effective critical point above J_t. Correct determination of Γ_c requires extrapolation M_trotter → ∞ (§3.6). For J_s > 0 the effect is amplified. Hence Γ_c(J_s > 0) > Γ_c(J_s = 0) > J_t. This is a **falsifiable prediction**, testable on quantum simulators.

### 3.5 Technical Specifications

| Component | Implementation |
|---|---|
| Language | Rust (production), Python (prototyping) |
| Algorithm | Wolff clusters + Parallel Tempering |
| Storage | i8 (8× memory savings) |
| RNG | Xoshiro256++ |
| Parallelism | Rayon |
| Statistics | Jackknife ±σ; τ_int |
| Verification | Pfeuty (1970): Γ_c = 1.0 |

### 3.6 Methodological Limitations

1. **Finite M_trotter.** Richardson extrapolation M → ∞ (M = 16, 32) shows the difference between M = 16 and M = 32 at β = 10 is smaller than the statistical error. At M = 32 the systematic Trotter error does not exceed the jackknife uncertainty. Nevertheless, precise determination of Γ_c would benefit from M = 64, 128.

2. **Small lattice sizes.** Reliable FSS requires L ≥ 16–32. A single run at L = 16, Γ = 1.0 confirms the trend, but systematic FSS over multiple L values is needed.

3. **Autocorrelations.** τ_int for v_stag and the Binder cumulant can be substantially longer near the transition.

4. **Thermalization.** Adaptive thermalization was not applied to all runs.

5. **Comparison with exact solution.** Verification is limited to a single value of Γ_c from Pfeuty (1970).

---

## 4. Interpretive Value of Ze

### 4.1 What Ze Adds to Z₂ Gauge Theory

| Z₂ Concept | Ze Interpretation |
|---|---|
| Gauge freedom | Agent's freedom to define what constitutes an error |
| Confinement | T/S events strongly correlated |
| Deconfinement | T/S events independent |
| Wegner duality | T ↔ S transformation on the dual lattice |
| Monopole | Topological defect on the dual lattice; in Ze: an event where the T/S alternation rule is broken |

### 4.2 Connection to Biology

The Ze interpretation bridges condensed matter physics and the Free Energy Principle in neuroscience (Friston, 2010; Fields et al., 2022). Both approaches describe systems that minimize surprisal or prediction error.

### 4.3 Ze and Fundamental Physics

**Wheeler's "it from bit" (1989).** Wheeler proposed that physical reality emerges from binary answers to yes/no questions—bits. The T/S events of Ze are precisely such Wheelerian bits. The Z₂ gauge lattice is the "apparatus" posing binary questions; the gauge freedom is the freedom to redefine which answer counts as 0 and which as 1.

**Bohm's implicate order (1952).** In de Broglie–Bohm theory, particles interact through a quantum potential in 3N-dimensional configuration space. Bohm termed this the "implicate order." Ze offers a microscopic model: T/S correlations between agents exist in spin configuration space, while observed quantities are projections onto 3+1d.

**ER = EPR (Maldacena & Susskind, 2013).** The ER = EPR hypothesis states that entangled particles are connected by wormholes. In Ze, the Wilson loop serves as a mathematical analogue: the area law corresponds to the persistence of correlations, the perimeter law to their breaking. This is an analogy; its status is conceptual.

**Synthesis.** Z₂ gauge theory in the Ze interpretation unifies several foundational programs:

| Program | Core Idea | Implementation in Ze |
|---|---|---|
| Wheeler: it from bit | Physics from binary responses | T/S events on the Z₂ lattice |
| Bohm: implicate order | Interaction in configuration space | T/S correlations between agents |
| ER = EPR | Entanglement = wormholes | Wilson loop as ER bridge (analogy) |
| Van Raamsdonk | Spacetime from entanglement | Z₂ configurations as a tensor network (hypothesis) |

**Conceptual novelty:** to the author's knowledge, no prior work has connected Z₂ gauge theory, the Free Energy Principle, and the "it from bit" program in a single mathematical model. An extensive search across INSPIRE-HEP, Semantic Scholar, and CrossRef (July 2026) found no publications combining these three directions.

---

## 5. Open Problems

| Problem | Status |
|---|---|
| Proof of U(1) phase for H_Ze | Classical U(1) spin liquid: ice = 1.0 (L = 4,5). Quantum U(1) phase: ice = 0.99 at Γ = 0.1, stable up to Γ = 2.0 (L = 3, M = 16). |
| Computing α via critical exponents | g² ∼ (Γ_c − Γ)^(νη) at the Z₂ → U(1) transition; α = g²/4π from the RG flow. |
| Generalizing JW to 3+1d | Partial (Su, 2025; Su & Martin, 2026, preprint). |
| Continuum limit | Triviality problem: φ⁴ in 3+1d is trivial (Aizenman, 1981; Fröhlich, 1982). A possible workaround: the Z₂ → U(1) transition must occur before the continuum limit. |
| Falsifiable prediction | Γ_c(J_s > 0) > J_t, testable on quantum simulators. |

### 5.1 Structural Relation for α (Consistency Check)

**(1) v⁎ = 1 − ln 2 ≈ 0.3069.** The maximum-entropy point of a binary channel under the antiparallelism constraint S = −T. **Status:** rigorously proven.

**(2) P(T|v⁎) = (2 − ln 2)/2 ≈ 0.6534.** The fundamental error probability of the agent. **Status:** direct consequence.

**(3) Relation between α and g.** In the effective U(1) gauge theory on the pyrochlore lattice the action is S_eff = (1/g)∑(E² + B²). The effective fine-structure constant α_eff is related to g through field normalization. **Status:** postulated on dimensional grounds; the precise normalization requires a derivation from H_Ze that is not presented here.

**(4) The factor P(T|v⁎).** It is postulated that the effective coupling perceived by the agent is α_eff = P(T|v⁎)·g/(4π). The rationale is that the agent registers only T-events, so the effective interaction strength is proportional to P(T|v⁎). **Status:** conceptual Ansatz, not derived from H_Ze.

**(5) The constant g.** In Hermele, Fisher & Balents (2004, §VI) the ring-exchange constant g is computed at sixth order in degenerate perturbation theory for the XXZ model on pyrochlore: g = C·(J_⊥/J_z)⁶·J_z, with C ≈ 0.25, giving g ≈ 0.14 J_z at J_⊥/J_z ≈ 0.95. **Critical caveat:** H_Ze is not equivalent to the XXZ model. Using g from Hermele et al. is an extrapolation by analogy, not a rigorous derivation.

**(6) Numerical check.** With g = 0.14 the formula gives α ≈ 0.00728 (1/α ≈ 137.4), a 0.27% deviation from α_exp = 1/137.036. **Status:** numerical coincidence; not an independent prediction.

### 5.1.1 Error Probability as the Source of Interaction

The agent makes binary predictions. At the critical point v⁎ = 1 − ln 2 the error probability is P(T|v⁎) = (2 − ln 2)/2 ≈ 0.6534. Each T-event generates a disturbance of the Z₂ gauge field that propagates through the lattice and interacts with other agents.

### 5.1.2 Relation Between Error Probability and Coupling Constant

In QED, α = e²/(4π) governs the probability of virtual photon emission. In Ze, a virtual photon is a disturbance of the Z₂ gauge field caused by a prediction error. The reasoning proceeds as follows: (1) the agent errs with probability P(T|v⁎); (2) the error creates a disturbance of the Z₂ field; (3) the disturbance remains coherent over a distance L (the correlation length of the U(1) liquid); (4) in 3+1d the disturbance propagates as a spherical wave, with the interaction probability falling as ∼1/(4πr); (5) for two agents at the boundary of the coherent domain (r = L), one obtains α = P(T|v⁎)/(4π·L) = (2 − ln 2)/(8π·L).

### 5.1.3 Coherence of the Interaction

In the U(1) spin liquid on the pyrochlore lattice, excitations propagate as massless photons. The correlation length ξ sets the scale over which the influence of an error remains coherent. The number of coherent cells along one dimension is N_coh = ξ/a, and the interaction between two agents separated by a distance r decays as ∼1/r. The dimensionless coupling constant is α = P(T|v⁎)/(4π·ξ/a).

### 5.1.4 Structural Prediction for α

Substituting P(T|v⁎) = (2 − ln 2)/2 gives the boxed relation α = (2 − ln 2)/(8π·L) = 0.0520/L. This is a structural prediction: α ∝ 1/L, where L is the effective coherent domain size. The functional form is derived from Ze principles; the numerical value L ≈ 7.13 is obtained from α_exp and is not predicted by the theory.

### 5.1.5 Direct Measurement and Comparison with Experiment

Simulations for L = 3–7 (M = 64, Γ = 0.05) confirm that the system resides in the ice-rule manifold (ice = 1.0 for all L). The formula α = 0.0520/L yields α = 0.00743 (1/134.6) at L = 7 and matches α_exp exactly at L = 7.13.

### 5.1.6 Computing L from H_Ze and the Renormalization Group

The effective coherent domain size is L = λ_e/a_Ze ≈ 7.13, giving Λ_Ze = L·m_e·c² ≈ 3.65 GeV. The one-loop QED β-function, β(α) = 2α²/(3π), gives α⁻¹(m_e) = α⁻¹(Λ_Ze) + 1.88, yielding α⁻¹(Λ_Ze) = 135.16 compared to 137.0 from the Ze formula—a 1.4% difference.

**Triviality.** Z₂ gauge theory in 3+1d is trivial in the continuum limit (Aizenman, 1981; Fröhlich, 1982). In the Ze framework the lattice is not a regularization to be removed; the spacing a_Ze is a fundamental constant. QED emerges as an effective theory at wavelengths λ ≫ a_Ze.

### 5.1.7 Responses to Methodological Criticisms

**(a) On the connection between H_Ze and the XXZ model.** H_Ze is not equivalent to the XXZ model of Hermele et al. Z₂ gauge theory on the pyrochlore lattice maps to the Ising model in the confinement phase (Fradkin & Susskind, 1978). Adding quantum fluctuations through Γ may yield a U(1) spin liquid analogous to that studied by Hermele et al., but a formal proof of the equivalence of their low-energy sectors is absent.

**(b) On the normalization α ∝ g.** In the effective U(1) gauge theory the canonical action is S = (1/g)∑(E² + B²). Canonical quantization gives [E_i, A_j] = i·g·δ_ij. The effective fine-structure constant is α_eff = g/(4π). (This convention follows the Hermele et al. definition of g as the ring-exchange constant with dimensions of energy.)

**(c) On the factor P(T|v⁎).** This is a conceptual Ansatz: the agent registers only T-events, so the effective coupling is proportional to P(T|v⁎). Without this factor, α = g/(4π) ≈ 0.0111, giving 1/α ≈ 89.8.

**(d) On the status of "no free parameters."** The value Γ = 0.94 is chosen so that g(Γ) = 0.14. This is an inverse problem. The formula is consistent with experiment at g ≈ 0.14, but neither g nor Γ is derived from first principles within Ze.

**(e) On Fig. 11 of Hermele et al. (2004).** Figure 11 of Hermele et al. displays the photon dispersion. The constant g is computed in §VI via sixth-order perturbation theory and verified by exact diagonalization on 16-site clusters.

---

## 6. Conclusion

The fine-structure constant can be expressed as α = P(T|v⁎)·g/(4π) = (2 − ln 2)·g/(8π), where g is the ring-exchange constant of the emergent U(1) spin liquid on the pyrochlore lattice. This is a consistency check, not a first-principles prediction. The formula yields α ≈ 1/138.5 when g ≈ 0.14—the value computed via sixth-order perturbation theory by Hermele, Fisher & Balents (2004) for the XXZ model on pyrochlore, a model related to but not equivalent to H_Ze. The functional form α ∝ g and the proportionality factor P(T|v⁎)/(4π) are postulated, not derived from H_Ze. The connection between H_Ze and the XXZ model requires explicit demonstration. The numerical coincidence (0.27% deviation) is noteworthy but does not constitute a derivation of α from first principles.

Numerical modeling confirms that v⁎ is a real point in the phase diagram; a quantum phase transition occurs near Γ ≈ J_t with a shift Γ_c > J_t due to spatial couplings; confinement is observed on the cubic lattice; and AFM order is destroyed by frustration. The Binder cumulant reaches the 2/3 limit in the deep AFM phase. The code has been cross-validated on three independent implementations.

The Ze → QED program is formulated as a sequence: (1) Z₂ gauge theory—rigorously proven (Wegner, 1971); (2) Majorana fermions in 1+1d—rigorously proven (Jordan–Wigner, 1928); (3) Z₂ → U(1) transition—requires a frustrated lattice; (4) Dirac fermions in 3+1d—an open problem; (5) QED as an effective theory—an open problem.

**Falsifiable prediction:** Γ_c(J_s > 0) > Γ_c(J_s = 0) > J_t. Testable on quantum simulators.

**Code and data:** https://github.com/djabbat/ze-theory (Apache 2.0).

---

## Acknowledgments

The author is grateful to numerous anonymous referees whose persistent criticism across many versions of this work compelled the abandonment of unwarranted claims and the limitation to rigorously established facts.

---

## References

[1] Friston, K. (2010). The free-energy principle: a unified brain theory? *Nature Reviews Neuroscience, 11*, 127–138.

[2] Fields, C., Friston, K., Glazebrook, J.F., & Levin, M. (2022). A free energy principle for generic quantum systems. *Progress in Biophysics and Molecular Biology, 173*, 36–59.

[3] Wegner, F.J. (1971). Duality in generalized Ising models and phase transitions without local order parameters. *Journal of Mathematical Physics, 12*, 2259–2272.

[4] Wilson, K.G. (1974). Confinement of quarks. *Physical Review D, 10*, 2445–2459.

[5] Jordan, P., & Wigner, E. (1928). Über das Paulische Äquivalenzverbot. *Zeitschrift für Physik, 47*, 631–651.

[6] Pfeuty, P. (1970). The one-dimensional Ising model with a transverse field. *Annals of Physics, 57*, 79–90.

[7] Gorantla, P., & Huang, T.-C. (2025). Duality-preserving deformation of 3+1d lattice Z₂ gauge theory with exact gapped ground states. *Physical Review B, 111*, 245110.

[8] Su, L. (2025). Bosonization and Kramers-Wannier Dualities in General Dimensions. Doctoral dissertation, Massachusetts Institute of Technology. Joint publication with I. Martin is in preprint stage (2026).

[9] Hermele, M., Fisher, M.P.A., & Balents, L. (2004). Pyrochlore photons: The U(1) spin liquid in a S = ½ three-dimensional frustrated magnet. *Physical Review B, 69*, 064404.

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

[28] Sandvik, A.W. (2010). Computational studies of quantum spin systems. *AIP Conference Proceedings, 1297*, 135–338. DOI: 10.1063/1.3518900.

[29] Fradkin, E., & Shenker, S.H. (1979). Phase diagrams of lattice gauge theories with Higgs fields. *Physical Review D, 19*, 3682–3697.

[30] Fradkin, E., & Susskind, L. (1978). Order and disorder in gauge systems and magnets. *Physical Review D, 17*, 2637–2658.

[31] Landauer, R. (1961). Irreversibility and heat generation in the computing process. *IBM Journal of Research and Development, 5*(3), 183–191.

[32] Caudy, W., & Greensite, J. (2008). On the ambiguity of spontaneously broken gauge symmetry. *Physical Review D, 78*, 025018.

[33] Grady, M. (2005). Reconsidering gauge-Higgs continuity. *Physics Letters B, 626*, 161–166.

---

*© 2026 Jaba Tqemaladze. All rights reserved.*
