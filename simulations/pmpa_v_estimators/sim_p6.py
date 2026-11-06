#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
P6 — category depletion and the quantum of birth (§7, protocol P6).

Tests the two readings of the information matrix M_C of a category.

  Archetype (working reading, Section 5.5):
      tau_c(0) = tau(M_C) is identical for every instance of the category,
      and initial budgets are quantized in units of tau0:
          tau(M_p) = n_p * tau0,  n_p in N.
      Observables: the proper-time budget of generation g — measured as the
      mean time-to-disintegration under identical environment — is stationary
      in g (slope ~ 0); no generation exhausts the category; and measured
      budgets across categories cluster at integer multiples of tau0.

  Stock (alternative reading):
      each birth depletes the matrix:  tau_g = tau0 * (1 - gamma*g),
      and past the point tau_g < tau_min births become dead — the category
      stops producing viable agents and its observable instances disappear.
      Observables: tau_g declines with g (slope significantly negative), and
      viable birth count collapses to zero at a finite generation.

The protocol is entirely in geometric reality: it needs only the statistics
of lifetimes against birth order, never a measurement of the pre-geometric
matrix itself.

Instance model: born with budget tau_g, spends proper time at a per-unit
rate R (lognormal, mean 1) — the standing cost of prediction error and
internal contacts — and disintegrates when the budget runs out:

    T_gi = tau_g / R_gi.

The mean lifetime of a generation is the measurable proxy for tau_g.

Part 2 (quantization): categories with primitive counts n_c = 1..K.  Under
the archetype+quantum reading the measured budgets cluster at multiples of
tau0; under a continuum reading they do not.

Run:
    python3 sim_p6.py --all
    python3 sim_p6.py --model stock --gamma 0.05 --jackknife
    python3 sim_p6.py --all --plot p6_birth.png

Only numpy is required (matplotlib optional, for --plot).
"""

import argparse
import math

import numpy as np

# ----------------------------------------------------------------------
# Spending-rate model (identical environment across generations)
# ----------------------------------------------------------------------
def spend_rate(n, rng, mu=1.0, sigma=0.35):
    """Per-unit proper-time spend rate R ~ lognormal(mean=mu)."""
    return rng.lognormal(mean=math.log(mu) - 0.5 * sigma * sigma, sigma=sigma,
                         size=n)


def lifetime(tau, R):
    """Time to disintegration: budget divided by spend rate."""
    return np.asarray(tau, dtype=float) / np.asarray(R, dtype=float)


# ----------------------------------------------------------------------
# Generation-sequence models
# ----------------------------------------------------------------------
def generation_budgets(model, g, tau0, gamma):
    """tau_g for generation g (0-indexed)."""
    if model == "archetype":
        return tau0
    if model == "stock":
        return tau0 * max(0.0, 1.0 - gamma * g)
    raise ValueError(model)


def simulate_category(model, tau0, tau_min, gamma, G, n_inst, rng):
    """Returns mean lifetime and viable count per generation.

    A birth with tau_g < tau_min is a dead birth: no viable agent, and the
    instance is absent from the observable sample.
    """
    means, viables, tau_g_list = [], [], []
    for g in range(G):
        tau_g = generation_budgets(model, g, tau0, gamma)
        tau_g_list.append(tau_g)
        if tau_g < tau_min:
            means.append(np.nan)
            viables.append(0)
            continue
        R = spend_rate(n_inst, rng)
        T = lifetime(tau_g, R)
        means.append(float(T.mean()))
        viables.append(int(n_inst))
    return np.array(means), np.array(viables), np.array(tau_g_list)


def slope_fit(g, means):
    """Least-squares slope of mean lifetime vs generation (viable gens)."""
    m = ~np.isnan(means)
    gv = g[m]
    yv = means[m]
    if len(gv) < 2:
        return float("nan"), float("nan")
    A = np.column_stack([np.ones_like(gv), gv])
    c, *_ = np.linalg.lstsq(A, yv, rcond=None)
    return float(c[0]), float(c[1])


def slope_jackknife(g, means):
    """Block-jackknife SE on the slope (delete-one-generation)."""
    m = ~np.isnan(means)
    gv = g[m]
    yv = means[m]
    n = len(gv)
    if n < 3:
        return float("nan")
    _, slope_full = slope_fit(gv, yv)
    slopes = np.zeros(n)
    for k in range(n):
        gk = np.delete(gv, k)
        yk = np.delete(yv, k)
        _, slopes[k] = slope_fit(gk, yk)
    var = (n - 1.0) / n * np.sum((slopes - slopes.mean()) ** 2)
    return float(slopes.mean()), float(math.sqrt(var)), slope_full


def report_generations(model, tau0, tau_min, gamma, G, n_inst, seed, jackknife):
    rng = np.random.default_rng(seed)
    g = np.arange(G)
    means, viables, tau_g = simulate_category(
        model, tau0, tau_min, gamma, G, n_inst, rng)

    if jackknife:
        _, slope_se, slope = slope_jackknife(g, means)
    else:
        _, slope = slope_fit(g, means)
        slope_se = float("nan")

    extinct = G
    for i, v in enumerate(viables):
        if v == 0:
            extinct = i
            break

    # discrimination: slope significant relative to jackknife SE, or > 1e-3
    sig = False
    if not np.isnan(slope_se):
        sig = bool(abs(slope) > 3.0 * slope_se)
    elif not np.isnan(slope):
        sig = bool(abs(slope) > 1e-3)

    alive = ~np.isnan(means)
    n_alive = int(alive.sum())
    if model == "archetype":
        verdict = "ARCHETYPE (stationary)" if (not sig and n_alive == G) \
            else "NOT ARCHETYPE"
    else:
        decline = sig or (extinct < G)
        verdict = "STOCK (depleting)" if decline else "NOT STOCK"

    print(
        f"\n=== P6 generation test [{model}] "
        f"(G={G}, n_inst/gen={n_inst}, tau0={tau0}, tau_min={tau_min}, "
        f"gamma={gamma}{', jackknife' if jackknife else ''}) ===")
    for i in range(G):
        v = viables[i]
        mn = f"{means[i]:8.3f}" if not np.isnan(means[i]) else "   dead "
        print(f"  g={i:3d}  tau_g={tau_g[i]:6.3f}  viable={v:4d}  "
              f"mean_life={mn}")
    print(f"  slope={slope:+.4f} (se={slope_se:.4f})  extinct_at_g={extinct}"
          f"  -> {verdict}")
    return slope, slope_se, extinct


# ----------------------------------------------------------------------
# Part 2 — quantization of initial budgets across categories
# ----------------------------------------------------------------------
def quantization_check(tau0, K, n_inst, seed, continuum=False):
    """Measured budgets vs n*tau0 arithmetic-progression test.

    archetype: tau(M_c) = c * tau0  => measured budgets form an arithmetic
    progression (a common multiplier <1/R>), so measured/c is constant
    continuum: tau(M_c) drawn continuously  => no regular spacing.
    Returns the coefficient of variation of measured/c.
    """
    rng = np.random.default_rng(seed)
    measured = []
    truth = []
    for c in range(1, K + 1):
        if continuum:
            tau_c = tau0 * (0.5 + 0.5 * rng.random() * K)   # continuous draw
        else:
            tau_c = c * tau0
        truth.append(tau_c)
        R = spend_rate(n_inst, rng)
        measured.append(float(lifetime(tau_c, R).mean()))
    measured = np.array(measured)
    truth = np.array(truth)
    ratio = measured / np.arange(1, K + 1)            # measured / c
    cv = float(ratio.std() / ratio.mean())
    print(
        f"\n=== P6 quantization [{'continuum' if continuum else 'archetype'}] "
        f"(K={K} categories, n_inst={n_inst}, tau0={tau0}) ===")
    for c in range(1, K + 1):
        print(f"  c={c:2d}  truth={truth[c-1]:6.3f}  measured={measured[c-1]:8.3f} "
              f"  measured/c={measured[c-1]/c:6.3f}")
    print(f"  CV(measured/c) = {cv:.4f} "
          f"(~0 => regular n*tau0 spacing; >0.1 => continuum)")
    return cv


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--all", action="store_true", help="both models + quantization")
    ap.add_argument("--model", choices=["archetype", "stock"])
    ap.add_argument("--G", type=int, default=12, help="generations")
    ap.add_argument("--n", type=int, default=300, help="instances per generation")
    ap.add_argument("--tau0", type=float, default=1.0)
    ap.add_argument("--tau_min", type=float, default=0.3,
                    help="viable-birth threshold (alpha*eps0)")
    ap.add_argument("--gamma", type=float, default=0.08, help="stock depletion rate")
    ap.add_argument("--K", type=int, default=8, help="categories for quantization")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--jackknife", action="store_true",
                    help="slope SE via delete-one-generation jackknife")
    ap.add_argument("--plot", metavar="OUT.png", help="save figure")
    a = ap.parse_args()

    if a.all:
        report_generations("archetype", a.tau0, a.tau_min, a.gamma, a.G, a.n,
                           a.seed, a.jackknife)
        report_generations("stock", a.tau0, a.tau_min, a.gamma, a.G, a.n,
                           a.seed, a.jackknife)
        rms_a = quantization_check(a.tau0, a.K, a.n, a.seed, continuum=False)
        rms_c = quantization_check(a.tau0, a.K, a.n, a.seed, continuum=True)
        _ = (rms_a, rms_c)
    elif a.model:
        report_generations(a.model, a.tau0, a.tau_min, a.gamma, a.G, a.n,
                           a.seed, a.jackknife)
    else:
        ap.print_help()
        return

    if a.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except Exception as e:  # pragma: no cover
            print(f"[plot skipped: {e}]")
            return
        rng = np.random.default_rng(a.seed)
        g = np.arange(a.G)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))
        for model, col in (("archetype", "tab:blue"), ("stock", "tab:red")):
            means, viables, tau_g = simulate_category(
                model, a.tau0, a.tau_min, a.gamma, a.G, a.n, rng)
            m = ~np.isnan(means)
            ax1.plot(g[m], means[m], "o-", color=col, ms=5, label=model)
        ax1.set_xlabel("generation g")
        ax1.set_ylabel("mean time-to-disintegration (proxy for tau_g)")
        ax1.set_title("P6: lifetime vs generation (archetype flat, stock declines)")
        ax1.legend()
        ax1.grid(alpha=0.3)

        for cont, col, lab in ((False, "tab:blue", "archetype (n*tau0)"),
                               (True, "tab:red", "continuum")):
            rng2 = np.random.default_rng(a.seed)
            meas = []
            for c in range(1, a.K + 1):
                if cont:
                    tau_c = a.tau0 * (0.5 + 0.5 * rng2.random() * a.K)
                else:
                    tau_c = c * a.tau0
                R = spend_rate(a.n, rng2)
                meas.append(lifetime(tau_c, R).mean())
            ax2.scatter(np.arange(1, a.K + 1), np.array(meas) / np.arange(1, a.K + 1),
                        color=col, s=42, label=lab)
        ax2.set_xlabel("category (primitive count c)")
        ax2.set_ylabel("measured budget / c")
        ax2.set_title("P6: quantization — measured/c flat for archetype")
        ax2.legend()
        ax2.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(a.plot, dpi=150)
        print(f"[plot saved: {a.plot}]")


if __name__ == "__main__":
    main()
