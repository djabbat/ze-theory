# Response to Peer Review — Ze Theory

**Author: Jaba Tqemaladze, MD | Date: 2026-07-04**

---

## Summary of Changes Made

The reviewer raised several valid criticisms. The following changes have been made:

1. **Abstract:** Clarified that g is computed via perturbation theory (not measured non-perturbatively); changed «no parameters fitted» to «consistency check».
2. **§5.1.0:** Completely rewritten to honestly characterize each step's logical status.
3. **§5.1.7 (NEW):** Direct responses to methodological criticisms.
4. **Conclusion:** Updated to reflect honest status.
5. **Reference [28] (Sandvik):** Fixed incorrect DOI.

---

## Point-by-Point Response

### Reviewer Point 1: Hermele Fig. 11 reference is incorrect

**Reviewer:** «g(Γ) is measured non-perturbatively through exact diagonalization (Hermele, Fisher & Balents, 2004, Fig. 11) — не соответствует содержанию статьи.»

**Response:** **ACCEPTED.** The reviewer is correct. Hermele et al. (2004) compute g via 6th-order degenerate perturbation theory, verified against exact diagonalization on 16-site clusters. The language «measured non-perturbatively» was wrong. Changes made:
- Abstract: replaced with «computed via 6th-order degenerate perturbation theory... verified against exact diagonalization»
- §5.1.0: added explicit description of the perturbative nature of g
- §5.1.7(e): acknowledged the error and corrected «Fig. 11» to «§VI»
- Removed all references to «Fig. 11» as a source of g measurement

### Reviewer Point 2: H_Ze ≠ XXZ model

**Reviewer:** «Связь между Z₂-калибровочной теорией и XXZ-моделью на пирохлорной решётке не продемонстрирована.»

**Response:** **ACCEPTED.** H_Ze is a Z₂ gauge theory with transverse field; the XXZ model on pyrochlore is a different Hamiltonian. Changes made:
- §5.1.0: added «КРИТИЧЕСКОЕ ЗАМЕЧАНИЕ: H_Ze НЕ эквивалентен XXZ-модели»
- §5.1.7(a): detailed the gap and suggested future work (duality mapping)
- Abstract: added «the connection between H_Ze and the XXZ model requires explicit demonstration which is not yet provided»

### Reviewer Point 3: 1/(4π) factor and α ∝ g

**Reviewer:** «Происхождение множителя 1/(4π) не обосновано. В стандартной теории α ∝ g², а не α ∝ g.»

**Response:** **PARTIALLY ACCEPTED.** 
- The 1/(4π) factor is a dimensional analysis Ansatz — this is now explicitly stated.
- Regarding α ∝ g vs α ∝ g²: In the effective U(1) gauge theory on pyrochlore, the canonical action is S = (1/g) Σ (E²+B²). After canonical quantization, [E_i, A_j] = i·g·δ_ij, and the effective fine-structure constant is α_eff = g/(2π) (see Hermele et al., Eqs. 6.8-6.10). This is α ∝ g, NOT α ∝ g². The reviewer's claim about α ∝ g² applies to standard QED where the coupling e² = 4πα appears in the covariant derivative. The emergent U(1) gauge theory on pyrochlore has different normalization.
- Changes: §5.1.7(b) added with explicit reference to Hermele Eqs. 6.8-6.10.

### Reviewer Point 4: P(T|v*) as multiplicative factor

**Reviewer:** «Это ad hoc введение вероятности ошибки как эффективного коэффициента ослабления связи.»

**Response:** **ACCEPTED.** This is the weakest link in the derivation. The factor P(T|v*) is a conceptual Ansatz, not a mathematical consequence of H_Ze. Changes made:
- §5.1.0(4): explicitly marked as «концептуальный Ansatz, не выведенный из H_Ze»
- §5.1.7(c): acknowledged that without this factor, α = g/(4π) = 0.0111 → 1/α = 89.8
- Conclusion: added «the proportionality factor P(T|v*)/(4π) [is] postulated based on dimensional analysis and conceptual Ansatz»

### Reviewer Point 5: Hidden fitting through L

**Reviewer:** «Если L = 7.13 получено подгонкой, то всё вычисление содержит подгонку.»

**Response:** **ACCEPTED.** In the original manuscript, L = 7.13 was indeed fitted to α_exp. In the revised version, the derivation has been restructured to use g directly (α = P(T|v*)·g/(4π)) rather than through L. However, Γ = 0.94 is still chosen to make g(Γ) = 0.14 (inverse problem). Changes made:
- Abstract: «consistency check, not a first-principles prediction» replaces «no parameters fitted»
- §5.1.0(6): explicitly states «Γ выбрано для получения g = 0.14»
- §5.1.7(d): honest characterization of the logical status

### Reviewer Point 6: Triviality problem

**Reviewer:** «Отказ от континуального предела не является решением проблемы тривиальности.»

**Response:** **PARTIALLY ACCEPTED.** The reviewer is correct that simply declaring the lattice fundamental does not solve the mathematical triviality problem. However, in an effective field theory framework, a finite cutoff is physically acceptable (e.g., the Standard Model is widely believed to have a finite cutoff at the Planck scale). The Z₂ lattice provides such a cutoff. Changes:
- This remains in the manuscript as stated, but the discussion has been framed as «effective theory» rather than «solution to triviality.»

### Reviewer Point 7: No falsifiable predictions

**Reviewer:** «Единственное предсказание... не специфицировано численно.»

**Response:** **PARTIALLY ACCEPTED.** The manuscript does contain one falsifiable prediction: Γ_c(J_s > 0) > Γ_c(J_s = 0) > J_t. However, the reviewer is correct that this is qualitative. Numerical specificity would strengthen it. Suggested improvement for future work: compute Γ_c(J_s = 0.1, 0.2, 0.3) with error bars from FSS.

---

## Updated Status of Key Claims

| Claim | Original Status | Reviewer Assessment | Revised Status |
|-------|----------------|---------------------|----------------|
| v* = 1−ln2 is max-entropy | Proven | Not disputed | ✅ Proven |
| H_Ze = Z₂ gauge theory | Proven | Not disputed | ✅ Proven |
| α = P(T\|v*)·g/(4π) | Derived | Ad hoc Ansatz | ⚠️ Postulated |
| g = 0.14 from Hermele | Measured via ED | Computed via PT | ⚠️ PT + ED verification |
| H_Ze ↔ XXZ on pyrochlore | Implicitly assumed | Not demonstrated | ❌ Open problem |
| No parameter fitting | Claimed | Hidden fitting | ⚠️ Consistency check |
| Triviality solved | Claimed | Not solved | ⚠️ Effective theory framing |

---

## Target Journal Recommendations

Given the honest reassessment:

**Primary target (with revisions):** *Foundations of Physics* (IF ~2.5)
- Welcomes interpretive frameworks
- Publishes conceptual work at the physics-philosophy interface
- The Ze framework fits the journal's scope

**Alternative:** *European Physical Journal B* (IF ~1.6)
- Condensed matter + interdisciplinary
- Numerical methodology is a strength here

**Not recommended:** *Physical Review Letters*, *Nature Physics*, *Physical Review X*
- Require rigorous derivations and/or experimental predictions
- The Ansatz-based nature of the α formula does not meet their standards

---

## Remaining Open Problems (for future work)

1. **Rigorous connection H_Ze → XXZ on pyrochlore:** Duality mapping or effective low-energy theory derivation.
2. **Independent determination of Γ*:** Physical principle (not inverse problem) that fixes Γ.
3. **Derivation of P(T|v*) factor:** Why does agent error probability multiply the gauge coupling?
4. **SSE on L≥8 pyrochlore:** Direct non-perturbative measurement of g in H_Ze (not XXZ model).
5. **Numerical specificity of Γ_c(J_s) prediction:** Quantitative values with error bars.
