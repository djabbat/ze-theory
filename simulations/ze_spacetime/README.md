# Ze System — Spacetime Emergence (L = 4..12)

Code for: **J. Tqemaladze**, *Information-Theoretic Emergence of Spacetime and
the Arrow of Time* (Table 1, Section 4).

A discrete information-theoretic toy model in which a quadratic form with
Lorentzian signature `Q = N_T² − N_S²` arises from counting two types of
binary events:

- **T-events (stasis)** — `x_{k+1} == x_k`
- **S-events (switch)** — `x_{k+1} != x_k`

## Files

| File | Purpose |
|------|---------|
| `ze_spacetime_enum.py` | Exhaustive enumeration of all `2^L` sequences, L = 4..12 |
| `table1.csv` | Machine-readable Table 1 (timelike/spacelike/null counts, entropy) |

## Usage

```bash
python3 ze_spacetime_enum.py              # print Table 1 (markdown)
python3 ze_spacetime_enum.py --self-test  # verify vs closed forms
python3 ze_spacetime_enum.py -L 4 12 -o table1.csv
```

## Results (verified)

- **Odd N** (= even L): timelike = spacelike = 2^N, null = 0 → exact 50/50 split.
- **Even N** (= odd L): null = 2·C(N, N/2), timelike = spacelike = 2^N − C(N, N/2).
- `S_max = log C(N, floor(N/2))` grows approximately linearly with N.

`--self-test` reproduces Table 1 exactly from these closed forms.

## Environment

Python 3.x, standard library only (`itertools`, `math`, `collections`).
Platform-independent. MIT License.
