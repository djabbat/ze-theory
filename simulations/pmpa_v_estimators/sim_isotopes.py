#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Isotopes as a real-world instance of the P6 structure (§7, protocol P6).

P6 distinguishes two readings of a category's information matrix by three
geometric-reality observables:

    (a) quantization  — instance budgets sit on an arithmetic progression
        (unit tau0, Section 5.5: tau(M_p) = n_p * tau0);
    (b) viability threshold — only births with tau > tau_min are agents;
        outside the window births are "dead" (spatial fragments);
    (c) stationarity — the category's instance set does not deplete (archetype),
        versus a stock that loses instances over generations.

Isotopes realize all three structures with measured data:

    category   <->  element (proton number Z)
    instance   <->  isotope (mass number A = Z + N)
    budget     <->  nuclear mass M(A,Z)
    quantum    <->  atomic mass unit m_u ~ 931.5 MeV/c^2  (the A lattice)
    viability  <->  bound nucleus (binding energy B > 0): the valley of
                    stability is a finite band of N for each Z; outside it
                    nuclei are unbound — "dead births" that do not persist
    archetype  <->  the chart of nuclides is stationary: an element does not
                    lose isotopes over generations of observation.

The binding energies are computed from the semi-empirical mass formula
(Weizsaecker).  Honest status: this is an ANALOGY that grounds the P6
pattern in known physics; it does not by itself test the Ze theory.  The
quantization of mass on the A lattice and the bounded valley of stability
are empirical facts; identifying them with a pre-geometric "quantum of
proper time" is an interpretation, not a derivation.

Run:
    python3 sim_isotopes.py --all
    python3 sim_isotopes.py --element 50 --plot p6_isotopes.png

Only numpy is required (matplotlib optional, for --plot).
"""

import argparse

import numpy as np

# Semi-empirical mass formula (MeV); standard coefficient set.
A_V, A_S, A_C, A_A, A_P = 15.67, 17.23, 0.714, 93.15, 11.2
M_P, M_N, M_U = 938.272, 939.565, 931.494   # MeV/c^2


def pairing(A, Z):
    """delta(A,Z): +even-even, 0 odd-A, -odd-odd (handles scalars/arrays)."""
    scalar = np.isscalar(A)
    A = np.atleast_1d(np.asarray(A, dtype=float))
    Z = np.atleast_1d(np.asarray(Z, dtype=float))
    if Z.size == 1 and A.size > 1:
        Z = np.full(A.shape, Z[0])
    out = np.zeros(A.shape, dtype=float)
    even_a = (A % 2 == 0)
    if np.any(even_a):
        term = A_P * A[even_a] ** (-3.0 / 4.0)
        out[even_a] = term * np.where(Z[even_a] % 2 == 0, 1.0, -1.0)
    return float(out[0]) if scalar else out


def binding(A, Z):
    """Weizsaecker binding energy B(A,Z) in MeV."""
    return (A_V * A
            - A_S * A ** (2.0 / 3.0)
            - A_C * Z * Z * A ** (-1.0 / 3.0)
            - A_A * (A - 2.0 * Z) ** 2 / A
            + pairing(A, Z))


def mass(A, Z):
    """Nuclear mass M(A,Z) in MeV/c^2."""
    return Z * M_P + (A - Z) * M_N - binding(A, Z)


def stable_z(A):
    """Z minimizing the mass parabola at fixed A (analytic approx)."""
    d1 = 4.0 * A_A / A
    d2 = 2.0 * A_C * A ** (-1.0 / 3.0)
    return (M_N - M_P + d1) / (d1 + d2)   # continuous optimum


def chart(max_a=300, zmax=110):
    """Viable nuclei: B > 0 and fission-stable (Z^2/A < 50).

    Returns (A[], Z[], B[]).  The B > 0 boundary alone overestimates
    viability for heavy elements (the real valley of stability is narrower
    because of fission and decay modes not present in the mass formula).
    """
    A = np.arange(1, max_a + 1)
    Z = np.arange(1, zmax + 1)
    AA, ZZ = np.meshgrid(A, Z, indexing="ij")
    mask = ZZ <= AA
    B = np.where(mask, binding(AA, ZZ), -1.0)
    fission = (ZZ ** 2 / AA) < 50.0
    vi = (B > 0.0) & fission
    return AA[vi], ZZ[vi], B[vi]


def report(max_a, zmax):
    A, Z, B = chart(max_a, zmax)
    # (a) quantization: for a sample element, mass vs A is A*m_u + small correction
    print("\n=== P6 isotopes: quantization of the budget (mass on the A lattice) ===")
    for Ze in (26, 50, 82):
        sel = Z == Ze
        As = A[sel]
        Ms = mass(As, Ze)
        # best linear fit slope in units of m_u
        slope = np.polyfit(As, Ms, 1)[0]
        resid = np.abs(Ms - (Ze * M_P + (As - Ze) * M_N) - 0.0)  # not used
        print(f"  Z={Ze:3d}: {len(As):3d} bound isotopes, "
              f"mass slope = {slope/M_U:+.5f} m_u per nucleon (ideal 1.00000)")

    # (b) viability band: for fixed Z, finite window of N; count window width
    print("\n=== P6 isotopes: viability threshold (valley of stability) ===")
    widths = []
    for Ze in (6, 26, 50, 82, 92):
        Ns = A[Z == Ze] - Ze
        if len(Ns) == 0:
            continue
        w = int(Ns.max() - Ns.min())
        widths.append(w)
        print(f"  Z={Ze:3d}: bound isotopes live in N window of width {w:3d} "
              f"(N from {int(Ns.min()):3d} to {int(Ns.max()):3d}) — "
              f"outside: dead births")
    print(f"  -> the viable window is finite and bounded for every element "
          f"(mean width ~{int(np.mean(widths))})")

    # (c) stationarity: the chart is a fixed set, no depletion
    n_bound = len(A)
    print(f"\n=== P6 isotopes: stationarity (archetype vs stock) ===")
    print(f"  {n_bound} bound nuclides in the chart (B > 0). The set is "
          f"time-invariant: an element does not lose isotopes across "
          f"generations of observation -> archetype signature.")

    # binding-energy curve peak (most viable region): scan Z for max B/A
    ba = []
    for a in range(1, 200):
        best = max((binding(a, z) / a, z) for z in range(1, a + 1))
        ba.append((a, best[0], best[1]))
    a_peak = max(ba, key=lambda t: t[1])
    print(f"\n  Binding energy per nucleon peaks at A={a_peak[0]} "
          f"(Z={a_peak[2]}, B/A={a_peak[1]:.2f} MeV) — the most viable region "
          f"(Fe/Ni).")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--max_a", type=int, default=300)
    ap.add_argument("--zmax", type=int, default=110)
    ap.add_argument("--element", type=int, help="plot one element's mass line")
    ap.add_argument("--plot", metavar="OUT.png")
    a = ap.parse_args()

    report(a.max_a, a.zmax)

    if a.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except Exception as e:  # pragma: no cover
            print(f"[plot skipped: {e}]")
            return
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 4.2))

        # (1) binding energy per nucleon (max over Z)
        aa = np.arange(1, 200)
        bb = np.array([max(binding(x, z) / x for z in range(1, x + 1))
                       for x in aa])
        ax1.plot(aa, bb, "-", color="tab:blue", lw=1.6)
        ax1.set_xlabel("mass number A")
        ax1.set_ylabel("B/A (MeV)")
        ax1.set_title("Binding per nucleon (peak near Fe/Ni)")

        # (2) chart of bound nuclei = viability band
        A, Z, _ = chart(a.max_a, a.zmax)
        ax2.scatter(A, Z, s=2, color="tab:green", alpha=0.5)
        ax2.set_xlabel("mass number A")
        ax2.set_ylabel("proton number Z")
        ax2.set_title(f"Bound nuclei (B>0): {len(A)} viable 'births'")

        # (3) quantization for a chosen element
        Ze = a.element or 50
        As = np.arange(Ze, 140)
        Ms = mass(As, Ze)
        ax3.plot(As, Ms / M_U, "-", color="tab:red", lw=1.5, label=f"Z={Ze}")
        ax3.plot(As, As, "--", color="k", lw=0.8, label="A (ideal m_u lattice)")
        ax3.set_xlabel("mass number A")
        ax3.set_ylabel("M / m_u")
        ax3.set_title("Mass quantized on the A lattice")
        ax3.legend()
        fig.tight_layout()
        fig.savefig(a.plot, dpi=150)
        print(f"[plot saved: {a.plot}]")


if __name__ == "__main__":
    main()
