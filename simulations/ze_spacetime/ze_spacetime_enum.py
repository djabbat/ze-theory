#!/usr/bin/env python3
"""
ze_spacetime_enum.py — Exhaustive enumeration of the Ze system state space.

Accompanies: J. Tqemaladze, "Information-Theoretic Emergence of Spacetime
and the Arrow of Time", Table 1 (Section 4).

Definitions (paper, Section 2):
  * T-event (stasis):   x_{k+1} == x_k
  * S-event (switch):   x_{k+1} != x_k
  * N = L - 1 transitions,  N_T + N_S = N
  * tau = sqrt(N_T^2 - N_S^2)                  (Ze proper time)
  * Z   = N_S / N_T                            (Ze impedance)
  * v   = N_S / N                              (Ze velocity)
  * S   = log C(N, N_S)                        (Boltzmann entropy; the
          constant log 2 from the free initial bit is omitted — it cancels
          in ratios and does not affect the maximum)

Classification (paper, Definition 6):
  * timelike:  N_T > N_S
  * spacelike: N_S > N_T
  * null:      N_T == N_S

Usage:
  python3 ze_spacetime_enum.py                 # print Table 1 (L = 4..12)
  python3 ze_spacetime_enum.py -L 4 12 -o table1.csv
  python3 ze_spacetime_enum.py --self-test     # verify against closed forms

Closed forms (used by --self-test):
  * odd  N: timelike = spacelike = 2^N, null = 0
  * even N: null = 2 * C(N, N/2),
            timelike = spacelike = 2^N - C(N, N/2)
  * S_max = log C(N, floor(N/2))

License: MIT. Standard library only (Python 3.x).
"""

import argparse
import itertools
import math
import sys
from collections import Counter


# ── Core enumeration ────────────────────────────────────────────────────

def analyze(L: int):
    """Exhaustively enumerate all 2^L binary sequences of length L.

    Returns dict with totals, timelike/spacelike/null counts, the maximal
    Boltzmann entropy, and the (N_T, N_S) histogram.
    """
    N = L - 1
    counts = {"timelike": 0, "spacelike": 0, "null": 0}
    ns_hist = Counter()

    for bits in itertools.product((0, 1), repeat=L):
        n_t = sum(1 for k in range(N) if bits[k + 1] == bits[k])
        n_s = N - n_t
        ns_hist[n_s] += 1
        if n_t > n_s:
            counts["timelike"] += 1
        elif n_s > n_t:
            counts["spacelike"] += 1
        else:
            counts["null"] += 1

    s_max = max(math.log(math.comb(N, ns)) for ns in ns_hist)
    return {
        "L": L, "N": N,
        "total": 2 ** L,
        **counts,
        "frac_timelike": counts["timelike"] / (2 ** L),
        "s_max": s_max,
        "hist": ns_hist,
    }


# ── Table rendering ─────────────────────────────────────────────────────

def render_table(rows, sep="|"):
    """Render the rows as the markdown Table 1 of the paper."""
    header = f"{sep} L (bits) {sep} N = L-1 {sep} Total {sep} Timelike {sep} Spacelike {sep} Null {sep} Frac. timelike {sep} Max S {sep}"
    rule = f"{sep}---{sep}---{sep}---{sep}---{sep}---{sep}---{sep}---{sep}---{sep}"
    lines = [header, rule]
    for r in rows:
        lines.append(
            f"{sep} {r['L']} {sep} {r['N']} {sep} {r['total']} {sep} "
            f"{r['timelike']} {sep} {r['spacelike']} {sep} {r['null']} {sep} "
            f"{r['frac_timelike']:.3f} {sep} {r['s_max']:.3f} {sep}"
        )
    return "\n".join(lines)


def write_csv(rows, path):
    """Write a machine-readable table (N_T, N_S, tau, Z, v, S per state)."""
    with open(path, "w", encoding="utf-8") as f:
        f.write("L,N,total,timelike,spacelike,null,frac_timelike,s_max\n")
        for r in rows:
            f.write(
                f"{r['L']},{r['N']},{r['total']},{r['timelike']},"
                f"{r['spacelike']},{r['null']},{r['frac_timelike']:.9f},"
                f"{r['s_max']:.9f}\n"
            )


# ── Self-test: closed forms ─────────────────────────────────────────────

def closed_form(L: int):
    """Analytic predictions for the counts (used for verification)."""
    N = L - 1
    if N % 2 == 1:                      # odd N: no null states, exact 50/50
        return 2 ** N, 2 ** N, 0
    null = 2 * math.comb(N, N // 2)     # even N: null = 2*C(N, N/2)
    rest = 2 ** L - null
    return rest // 2, rest // 2, null


def self_test(L_min=4, L_max=12):
    """Verify enumeration against the closed-form expressions."""
    failures = 0
    for L in range(L_min, L_max + 1):
        r = analyze(L)
        tl, sp, nu = closed_form(L)
        s_max = max(
            math.log(math.comb(r["N"], ns))
            for ns in range(r["N"] + 1)
        )
        ok = (r["timelike"] == tl and r["spacelike"] == sp
              and r["null"] == nu and abs(r["s_max"] - s_max) < 1e-9)
        status = "ok" if ok else "FAIL"
        if not ok:
            failures += 1
        print(f"  L={L:>2}  N={r['N']:>2}  total={r['total']:>5}  "
              f"timelike={r['timelike']:>5}  spacelike={r['spacelike']:>5}  "
              f"null={r['null']:>5}  frac={r['frac_timelike']:.6f}  "
              f"s_max={r['s_max']:.6f}  [{status}]")
    print()
    if failures:
        print(f"❌ {failures} L-value(s) FAILED")
        return 1
    print("✅ All closed forms match the enumeration exactly.")
    return 0


# ── CLI ─────────────────────────────────────────────────────────────────

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-L", nargs=2, type=int, default=[4, 12],
                    metavar=("MIN", "MAX"), help="bit-length range (default 4 12)")
    ap.add_argument("-o", "--output", metavar="CSV", help="write table to CSV")
    ap.add_argument("--self-test", action="store_true",
                    help="verify enumeration against closed forms")
    args = ap.parse_args(argv)

    if args.self_test:
        return self_test(args.L[0], args.L[1])

    rows = [analyze(L) for L in range(args.L[0], args.L[1] + 1)]
    print(render_table(rows))
    if args.output:
        write_csv(rows, args.output)
        print(f"\nCSV written to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
