#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PMPA overdetermination test — measure the Ze velocity v three ways (§7.6).

Implements the operational protocol of "The Ze Principle of Mutual Predictive
Adaptation" (§7.6): the Ze velocity v = K^2 (which-path knowledge) must be
recoverable from three independent observable channels of one two-slit run:

    v1  detector side : v1 = 1 - o_est^2,  o_est = measured marker-state
                        fidelity (calibrated path-conditional readouts)
    v2  screen  side  : v2 = 1 - (V_meas/V0)^2,  V_meas = fringe contrast
                        from the screen histogram
    v3  eraser  side  : v3 = 1 - (2 P+ - 1)^2,  P+ = fraction of idler
                        outcomes in the erased (|+>) basis (delayed-choice
                        coincidence statistics)

Two physical models are compared.

  Model A (coherent / Englert-saturating; standard QM):
        |Psi> = (|psi_A>|D_A> + |psi_B>|D_B>)/sqrt(2),  <D_A|D_B> = o.
        Fringe visibility V = V0 * o  =>  v = 1 - o^2,  Englert-saturated.

  Model B (naive materialization; literal reading of Postulate I):
        each event is "corpuscular" (which-path, prob v, lands classically)
        or "wave" (prob 1-v, lands coherently).

The overdetermination test is a *discriminator*: under Model A the three
estimators must coincide (v1 ~ v2 ~ v3 = 1-o^2); under Model B they must
deviate from the model's own materialization parameter v, flagging the model
as internally inconsistent.  This is the sense in which the protocol is
falsifiable rather than a tautology.

Statistical errors: block-jackknife (delete-one-block resampling).  The run
is split into R independent blocks; each estimator is recomputed on the pooled
data with one block removed, and the standard error follows the jackknife
variance formula  se^2 = (R-1)/R * sum_k (theta_k - mean)^2.

Run:
    python3 sim_pmpa.py --all                     # both grids, full report
    python3 sim_pmpa.py --all --jackknife --blocks 24
    python3 sim_pmpa.py --model A --o 0.5         # single Model-A point
    python3 sim_pmpa.py --model B --v 0.3
    python3 sim_pmpa.py --all --jackknife --plot pmpa_visibility.png

Only numpy is required.
"""

import argparse
import math

import numpy as np

# ----------------------------------------------------------------------
# Geometry (dimensionless).  lambda*L = 1  =>  x is in units of lambda*L.
# ----------------------------------------------------------------------
LAMBDA_L = 1.0
SLIT_WIDTH = 0.2          # a  (envelope: zeros at |x| = 1/a = 5)
SLIT_SEP = 0.5            # d  (fringes: spacing = lambda*L/d = 2)
XMAX = 4.9
NPIX = 201


def grid():
    return np.linspace(-XMAX, XMAX, NPIX)


def edges():
    """Histogram edges whose bin centers coincide with the grid points."""
    g = grid()
    h = g[1] - g[0]
    return np.concatenate([[g[0] - h / 2.0], (g[:-1] + g[1:]) / 2.0, [g[-1] + h / 2.0]])


def envelope(x):
    """Single-slit intensity envelope |psi(x)|^2 (Fraunhofer)."""
    u = math.pi * SLIT_WIDTH * x / LAMBDA_L
    out = np.ones_like(x)
    m = np.abs(u) > 1e-12
    out[m] = (np.sin(u[m]) / u[m]) ** 2
    return out


def delta_phase(x):
    """Path phase difference between the two slits."""
    return 2.0 * math.pi * SLIT_SEP * x / LAMBDA_L


def pdf_int(x):
    """Coherent (full-interference) distribution, normalized over the grid."""
    p = envelope(x) * (1.0 + np.cos(delta_phase(x)))
    return p / p.sum()


def pdf_mix(x):
    """Classical (which-path) distribution: single-slit envelope only."""
    e = envelope(x)
    return e / e.sum()


def pdf_anti(x):
    """Anti-fringe pattern (idler |-> subset)."""
    p = envelope(x) * (1.0 - np.cos(delta_phase(x)))
    return p / p.sum()


def pdf_modelA(x, o):
    """Model A: p(x) ~ env * (1 + o cos delta), normalized."""
    p = envelope(x) * (1.0 + o * np.cos(delta_phase(x)))
    return p / p.sum()


def pdf_modelB(x, v):
    """Model B: p(x) = v * mix + (1-v) * int, normalized."""
    return v * pdf_mix(x) + (1.0 - v) * pdf_int(x)


# ----------------------------------------------------------------------
# Measurement helpers
# ----------------------------------------------------------------------
def sample_pdf(pdf_arr, n, rng):
    idx = rng.choice(NPIX, size=n, p=pdf_arr)
    return grid()[idx]


def measure_visibility(hist):
    """Fit p(x) = env(x)*(c0 + c1*cos(delta)) on the nonzero-envelope region.

    Returns V = c1/c0 (fringe contrast).  V0 = 1 for this ideal setup.
    """
    x = grid()
    e = envelope(x)
    cd = np.cos(delta_phase(x))
    m = e > 1e-4 * e.max()
    A = np.column_stack([e[m], e[m] * cd[m]])
    c, *_ = np.linalg.lstsq(A, hist[m] / hist.sum(), rcond=None)
    return float(c[1] / c[0])


def o_est_from_arrays(rA, rB, nbins=120):
    """Marker-state fidelity from two path-conditional readout arrays.

    Bhattacharyya coefficient of the two empirical densities.
    """
    rA = np.asarray(rA, dtype=float)
    rB = np.asarray(rB, dtype=float)
    lo = min(rA.min(), rB.min())
    hi = max(rA.max(), rB.max())
    bins = np.linspace(lo, hi, nbins + 1)
    bw = bins[1] - bins[0]
    hA, _ = np.histogram(rA, bins=bins, density=True)
    hB, _ = np.histogram(rB, bins=bins, density=True)
    return float(bw * np.sum(np.sqrt(hA * hB)))


def estimate_from_pool(hist, n_plus, n_events, o_est):
    """(v1, v2, v3, V) from pooled run data."""
    V = measure_visibility(hist)
    v2 = 1.0 - V * V
    v3 = 1.0 - (2.0 * (n_plus / n_events) - 1.0) ** 2
    v1 = 1.0 - o_est ** 2
    return v1, v2, v3, V


# ----------------------------------------------------------------------
# Model A generators  (coherent / Englert-saturating)
# ----------------------------------------------------------------------
class ModelA:
    def __init__(self, o, seed):
        self.o = float(o)
        self.rng = np.random.default_rng(seed)
        self.sigma = 1.0
        if self.o <= 1e-9:
            self.mu = 50.0  # effectively orthogonal marker states
        else:
            self.mu = self.sigma * math.sqrt(-2.0 * math.log(self.o))

    def v_true(self):
        return 1.0 - self.o ** 2

    def run_block(self, n):
        """Returns (screen_histogram, n_plus, n_events)."""
        x = grid()
        pplus = (1.0 + self.o) / 2.0
        plus = self.rng.random(n) < pplus
        n_plus = int(plus.sum())
        x_plus = sample_pdf(pdf_int(x), n_plus, self.rng)
        x_minus = sample_pdf(pdf_anti(x), n - n_plus, self.rng)
        hist, _ = np.histogram(
            np.concatenate([x_plus, x_minus]), bins=edges()
        )
        return hist, n_plus, n

    def run(self, n_events):
        hist, n_plus, n = self.run_block(n_events)
        return float(n_plus) / n, hist

    def calibration_arrays(self, path, n):
        s = 1.0 if path == "A" else -1.0
        return self.rng.normal(s * self.mu, self.sigma, size=n)


# ----------------------------------------------------------------------
# Model B generators  (naive materialization of Postulate I)
# ----------------------------------------------------------------------
class ModelB:
    def __init__(self, v, seed):
        self.v = float(v)
        self.rng = np.random.default_rng(seed)
        self.sigma = 1.0
        self.mu = 8.0

    def v_true(self):
        return self.v

    def run_block(self, n):
        """Returns (screen_histogram, n_plus, n_events)."""
        x = grid()
        n_corp = int(self.rng.binomial(n, self.v))
        n_wave = n - n_corp
        if n_corp:
            plus_c = self.rng.random(n_corp) < 0.5
            x_corp = sample_pdf(pdf_mix(x), n_corp, self.rng)
        else:
            plus_c = np.array([], dtype=bool)
            x_corp = np.array([])
        if n_wave:
            plus_w = np.ones(n_wave, dtype=bool)
            x_wave = sample_pdf(pdf_int(x), n_wave, self.rng)
        else:
            plus_w = np.array([], dtype=bool)
            x_wave = np.array([])
        n_plus = int(np.concatenate([plus_c, plus_w]).sum())
        hist, _ = np.histogram(
            np.concatenate([x_corp, x_wave]), bins=edges()
        )
        return hist, n_plus, n

    def run(self, n_events):
        hist, n_plus, n = self.run_block(n_events)
        return float(n_plus) / n, hist

    def calibration_arrays(self, path, n):
        """Path-conditional readouts: v fraction separated, 1-v ambiguous."""
        s = 1.0 if path == "A" else -1.0
        n_sep = int(self.rng.binomial(n, self.v))
        sep = self.rng.normal(s * self.mu, self.sigma, size=n_sep)
        amb = self.rng.normal(0.0, self.sigma, size=n - n_sep)
        return np.concatenate([sep, amb])


# ----------------------------------------------------------------------
# Estimators: single run and block-jackknife
# ----------------------------------------------------------------------
def estimate_single(model, n_events, n_cal):
    hist, n_plus, n = model.run_block(n_events)
    o_est = o_est_from_arrays(
        model.calibration_arrays("A", n_cal),
        model.calibration_arrays("B", n_cal),
    )
    return np.array(estimate_from_pool(hist, n_plus, n, o_est))


def jackknife_estimate(model, n_blocks, events_per_block, n_cal_per_block):
    """Block-jackknife: full-sample estimate + jackknife mean and SE."""
    R = n_blocks
    hists, nplus, nev = [], [], []
    rA, rB = [], []
    for _ in range(R):
        h, np_, n_ = model.run_block(events_per_block)
        hists.append(h)
        nplus.append(np_)
        nev.append(n_)
        rA.append(model.calibration_arrays("A", n_cal_per_block))
        rB.append(model.calibration_arrays("B", n_cal_per_block))

    full_hist = np.sum(hists, axis=0)
    full_np = sum(nplus)
    full_n = sum(nev)
    full_o = o_est_from_arrays(np.concatenate(rA), np.concatenate(rB))
    theta_full = np.array(
        estimate_from_pool(full_hist, full_np, full_n, full_o)
    )

    theta_k = np.zeros((R, 4))
    for k in range(R):
        h = full_hist - hists[k]
        np_ = full_np - nplus[k]
        n_ = full_n - nev[k]
        rA_k = np.concatenate([rA[i] for i in range(R) if i != k])
        rB_k = np.concatenate([rB[i] for i in range(R) if i != k])
        o = o_est_from_arrays(rA_k, rB_k)
        theta_k[k] = estimate_from_pool(h, np_, n_, o)

    mean = theta_k.mean(axis=0)
    var = (R - 1.0) / R * np.sum((theta_k - mean) ** 2, axis=0)
    se = np.sqrt(var)
    return theta_full, mean, se


def consistency(mean, v_true, tol=0.05):
    """PASS if all three estimators sit within tol of v_true."""
    ok = bool(np.all(np.abs(mean - v_true) <= tol))
    spread = float(mean.max() - mean.min())
    return ok, spread


# ----------------------------------------------------------------------
# CLI / reports
# ----------------------------------------------------------------------
NAMES = ["v1", "v2", "v3", "V"]


def report_point(model, label, n_events, n_cal, tol, jackknife, blocks):
    if jackknife:
        theta_full, mean, se = jackknife_estimate(
            model, blocks, n_events // blocks, n_cal // blocks
        )
        ok, spread = consistency(mean[:3], model.v_true(), tol)
        s = "  ".join(
            f"{NAMES[i]}={mean[i]:.3f}±{se[i]:.3f}" for i in range(4)
        )
        print(
            f"{label:11s} v_true={model.v_true():6.3f} | {s} | "
            f"spread={spread:5.3f} -> {'PASS' if ok else 'FAIL (v_true mismatch)'}"
        )
        return theta_full, mean, se, ok
    theta = estimate_single(model, n_events, n_cal)
    ok, spread = consistency(theta[:3], model.v_true(), tol)
    s = "  ".join(f"{NAMES[i]}={theta[i]:6.3f}" for i in range(4))
    print(
        f"{label:11s} v_true={model.v_true():6.3f} | {s} | "
        f"spread={spread:5.3f} -> {'PASS' if ok else 'FAIL (v_true mismatch)'}"
    )
    return theta, None, None, ok


def run_grid(model_name, n_events, n_cal, tol, seed, jackknife, blocks):
    print(
        f"\n=== Model {model_name}: v1 vs v2 vs v3 "
        f"({'jackknife R=' + str(blocks) if jackknife else 'single'}; "
        f"N~{n_events}, tol={tol}) ==="
    )
    rows = []
    if model_name == "A":
        grid_pts = [0.0, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, 1.0]
        for o in grid_pts:
            rows.append(report_point(
                ModelA(o, seed=seed), f"A(o={o:.2f})", n_events, n_cal, tol,
                jackknife, blocks))
    else:
        grid_pts = [0.0, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, 0.99]
        for v in grid_pts:
            rows.append(report_point(
                ModelB(v, seed=seed), f"B(v={v:.2f})", n_events, n_cal, tol,
                jackknife, blocks))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--model", choices=["A", "B"])
    ap.add_argument("--o", type=float, help="Model A: marker overlap")
    ap.add_argument("--v", type=float, help="Model B: materialization rate")
    ap.add_argument("--all", action="store_true", help="run both grids")
    ap.add_argument("--N", type=int, default=200000, help="events per point")
    ap.add_argument("--ncal", type=int, default=40000, help="calibration reads")
    ap.add_argument("--tol", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--jackknife", action="store_true", help="block-jackknife SE")
    ap.add_argument("--blocks", type=int, default=24, help="jackknife blocks")
    ap.add_argument("--plot", metavar="OUT.png", help="visibility-law plot")
    a = ap.parse_args()

    if a.plot or a.all:
        run_grid("A", a.N, a.ncal, a.tol, a.seed, a.jackknife, a.blocks)
        run_grid("B", a.N, a.ncal, a.tol, a.seed, a.jackknife, a.blocks)
    elif a.model == "A":
        report_point(ModelA(a.o, seed=a.seed), f"A(o={a.o:.2f})", a.N,
                     a.ncal, a.tol, a.jackknife, a.blocks)
    elif a.model == "B":
        report_point(ModelB(a.v, seed=a.seed), f"B(v={a.v:.2f})", a.N,
                     a.ncal, a.tol, a.jackknife, a.blocks)
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
        vs = np.linspace(0, 1, 200)
        V_A = np.sqrt(1.0 - vs)                     # Englert-saturated
        V_B = 2.0 * (1.0 - vs) / (2.0 - vs)         # Model B materialization
        fig, ax = plt.subplots(figsize=(6.4, 4.6))
        ax.plot(vs, V_A, "-", lw=1.6, label="V0*sqrt(1-v)  (Model A, QM)")
        ax.plot(vs, V_B, "--", lw=1.6, label="2(1-v)/(2-v)  (Model B)")
        err = a.jackknife
        for o in [0.1, 0.25, 0.5, 0.75, 0.9]:
            m = ModelA(o, seed=a.seed)
            if err:
                tf, mean, se = jackknife_estimate(
                    m, a.blocks, a.N // a.blocks, a.ncal // a.blocks)
                vt, V, seV = m.v_true(), mean[3], se[3]
            else:
                theta = estimate_single(m, a.N, a.ncal)
                vt, V, seV = m.v_true(), theta[3], 0.0
            ax.errorbar(vt, V, yerr=seV, fmt="o", color="tab:blue", ms=5)
        for v in [0.1, 0.25, 0.5, 0.75, 0.9]:
            m = ModelB(v, seed=a.seed)
            if err:
                tf, mean, se = jackknife_estimate(
                    m, a.blocks, a.N // a.blocks, a.ncal // a.blocks)
                vt, V, seV = m.v_true(), mean[3], se[3]
            else:
                theta = estimate_single(m, a.N, a.ncal)
                vt, V, seV = m.v_true(), theta[3], 0.0
            ax.errorbar(vt, V, yerr=seV, fmt="s", color="tab:red", ms=5)
        ax.set_xlabel("Ze velocity v")
        ax.set_ylabel("measured fringe visibility V")
        ax.set_title("Visibility law: Model A (QM) vs Model B (materialization)")
        ax.legend()
        fig.tight_layout()
        fig.savefig(a.plot, dpi=150)
        print(f"[plot saved: {a.plot}]")


if __name__ == "__main__":
    main()
