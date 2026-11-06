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
diverge, flagging the model as internally inconsistent.  This is the sense
in which the protocol is falsifiable rather than a tautology.

Run:
    python3 sim_pmpa.py --all                 # both grids, full report
    python3 sim_pmpa.py --model A --o 0.5     # single Model-A point
    python3 sim_pmpa.py --model B --v 0.3     # single Model-B point
    python3 sim_pmpa.py --all --plot pmpa.png # grids + visibility-law plot

Only numpy is required.
"""

import argparse
import math
import sys

import numpy as np

# ----------------------------------------------------------------------
# Geometry (dimensionless).  lambda*L = 1  =>  x is in units of lambda*L.
# ----------------------------------------------------------------------
LAMBDA_L = 1.0
SLIT_WIDTH = 0.2          # a  (envelope: zeros at |x| = 1/a = 5)
SLIT_SEP = 0.5            # d  (fringes: spacing = 1/d ... = lambda*L/d)
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
def measure_visibility(hist):
    """Fit p(x) = env(x)*(c0 + c1*cos(delta)) on the nonzero-envelope region.

    Returns V = c1/c0 (fringe contrast).  V0 = 1 for this ideal setup.
    """
    x = grid()
    e = envelope(x)
    cd = np.cos(delta_phase(x))
    m = e > 1e-4 * e.max()
    # least squares: design matrix [e, e*cos]
    A = np.column_stack([e[m], e[m] * cd[m]])
    c, *_ = np.linalg.lstsq(A, hist[m] / hist.sum(), rcond=None)
    return c[1] / c[0]


def sample_pdf(pdf_arr, n, rng):
    idx = rng.choice(NPIX, size=n, p=pdf_arr)
    return grid()[idx]


def calibration_overlap(rng, gen_readout, n_cal=20000, nbins=120):
    """Measure the marker-state fidelity from calibration readouts.

    gen_readout(path) -> array of readout samples for a known path.
    Returns o_est = Bhattacharyya coefficient of the two empirical
    path-conditional readout distributions.
    """
    rA = np.asarray(gen_readout("A"))
    rB = np.asarray(gen_readout("B"))
    lo = min(rA.min(), rB.min())
    hi = max(rA.max(), rB.max())
    bins = np.linspace(lo, hi, nbins + 1)
    bw = bins[1] - bins[0]
    hA, _ = np.histogram(rA, bins=bins, density=True)
    hB, _ = np.histogram(rB, bins=bins, density=True)
    # Bhattacharyya coefficient of the two densities (needs the bin width)
    return float(bw * np.sum(np.sqrt(hA * hB)))


# ----------------------------------------------------------------------
# Model A generators  (coherent / Englert-saturating)
# ----------------------------------------------------------------------
class ModelA:
    def __init__(self, o, seed):
        self.o = float(o)
        self.rng = np.random.default_rng(seed)
        # marker readout separation consistent with overlap o:
        #   o = exp(-mu^2/(2 sigma^2))  =>  mu = sigma*sqrt(-2 ln o)
        self.sigma = 1.0
        if self.o <= 1e-9:
            self.mu = 50.0  # effectively orthogonal marker states
        else:
            self.mu = self.sigma * math.sqrt(-2.0 * math.log(self.o))

    def v_true(self):
        return 1.0 - self.o ** 2

    def run(self, n_events):
        """One run: returns (idler_plus_fraction, screen_histogram)."""
        x = grid()
        pplus = (1.0 + self.o) / 2.0
        plus = self.rng.random(n_events) < pplus
        n_plus = int(plus.sum())
        # screen positions for the two idler subsets
        x_plus = sample_pdf(pdf_int(x), n_plus, self.rng)
        x_minus = sample_pdf(pdf_anti(x), n_events - n_plus, self.rng)
        hist, _ = np.histogram(
            np.concatenate([x_plus, x_minus]), bins=edges()
        )
        return float(n_plus) / n_events, hist

    def calibration_readout(self, path, n=20000):
        s = 1.0 if path == "A" else -1.0
        return self.rng.normal(s * self.mu, self.sigma, size=n)


# ----------------------------------------------------------------------
# Model B generators  (naive materialization of Postulate I)
# ----------------------------------------------------------------------
class ModelB:
    def __init__(self, v, seed):
        self.v = float(v)
        self.rng = np.random.default_rng(seed)
        # corpuscular marker: large separation (o_c ~ 0); wave marker: none
        self.sigma = 1.0
        self.mu = 8.0

    def v_true(self):
        return self.v

    def run(self, n_events):
        """Returns (idler_plus_fraction, screen_histogram)."""
        x = grid()
        n_corp = int(self.rng.binomial(n_events, self.v))
        n_wave = n_events - n_corp
        if n_corp:
            # corpuscular: definite path -> 50/50 idler, classical landing
            path = self.rng.integers(0, 2, size=n_corp)
            plus_c = self.rng.random(n_corp) < 0.5
            x_corp = sample_pdf(pdf_mix(x), n_corp, self.rng)
        else:
            plus_c = np.array([], dtype=bool)
            x_corp = np.array([])
        if n_wave:
            # wave: coherent -> idler always in |+> subset, fringe landing
            plus_w = np.ones(n_wave, dtype=bool)
            x_wave = sample_pdf(pdf_int(x), n_wave, self.rng)
        else:
            plus_w = np.array([], dtype=bool)
            x_wave = np.array([])
        n_plus = int(np.concatenate([plus_c, plus_w]).sum())
        hist, _ = np.histogram(
            np.concatenate([x_corp, x_wave]), bins=edges()
        )
        return float(n_plus) / n_events, hist

    def calibration_readout(self, path, n=20000):
        """Path-conditional readouts: v fraction separated, 1-v ambiguous."""
        s = 1.0 if path == "A" else -1.0
        n_sep = int(self.rng.binomial(n, self.v))
        sep = self.rng.normal(s * self.mu, self.sigma, size=n_sep)
        amb = self.rng.normal(0.0, self.sigma, size=n - n_sep)
        return np.concatenate([sep, amb])


# ----------------------------------------------------------------------
# Estimators (only measured data is used; the true path is never touched)
# ----------------------------------------------------------------------
def estimate(model, n_events, n_cal=40000):
    """Return (v1, v2, v3, V_meas) from one run + calibration."""
    pplus, hist = model.run(n_events)
    v3 = 1.0 - (2.0 * pplus - 1.0) ** 2                      # eraser side
    V = measure_visibility(hist)
    v2 = 1.0 - V * V                                         # screen side (V0=1)
    o_est = calibration_overlap(model.rng, model.calibration_readout, n_cal)
    v1 = 1.0 - o_est ** 2                                    # detector side
    return v1, v2, v3, V


def consistency(v1, v2, v3, v_true, tol=0.05):
    """PASS if all three estimators sit within tol of v_true and of each other."""
    vals = np.array([v1, v2, v3])
    close_to_true = np.all(np.abs(vals - v_true) <= tol)
    spread = vals.max() - vals.min()
    return close_to_true and spread <= tol, float(spread)


# ----------------------------------------------------------------------
# CLI / reports
# ----------------------------------------------------------------------
def report_point(model, label, n_events, n_cal, tol):
    v1, v2, v3, V = estimate(model, n_events, n_cal)
    vt = model.v_true()
    ok, spread = consistency(v1, v2, v3, vt, tol)
    print(
        f"{label:9s} v_true={vt:6.3f} | v1={v1:6.3f} v2={v2:6.3f} "
        f"v3={v3:6.3f} | V={V:6.3f} | spread={spread:5.3f} -> "
        f"{'PASS' if ok else 'FAIL (v_true mismatch)'}"
    )
    return vt, v1, v2, v3, V, ok


def run_grid(model_name, n_events, n_cal, tol, rng_seed):
    print(f"\n=== Model {model_name}: v1 vs v2 vs v3 (N={n_events}, tol={tol}) ===")
    rows = []
    if model_name == "A":
        for o in [0.0, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, 1.0]:
            m = ModelA(o, seed=rng_seed)
            rows.append(report_point(m, f"A(o={o:.2f})", n_events, n_cal, tol))
    else:
        for v in [0.0, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, 0.99]:
            m = ModelB(v, seed=rng_seed)
            rows.append(report_point(m, f"B(v={v:.2f})", n_events, n_cal, tol))
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
    ap.add_argument("--plot", metavar="OUT.png", help="visibility-law plot")
    a = ap.parse_args()

    if a.plot or a.all:
        run_grid("A", a.N, a.ncal, a.tol, a.seed)
        run_grid("B", a.N, a.ncal, a.tol, a.seed)
    elif a.model == "A":
        m = ModelA(a.o, seed=a.seed)
        report_point(m, f"A(o={a.o:.2f})", a.N, a.ncal, a.tol)
    elif a.model == "B":
        m = ModelB(a.v, seed=a.seed)
        report_point(m, f"B(v={a.v:.2f})", a.N, a.ncal, a.tol)
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
        # visibility law: QM/Model A vs naive materialization/Model B
        vs = np.linspace(0, 1, 200)
        V_A = np.sqrt(1.0 - vs)                     # Englert-saturated V0*sqrt(1-v)
        V_B = 2.0 * (1.0 - vs) / (2.0 - vs)         # Model B materialization law
        fig, ax = plt.subplots(figsize=(6.4, 4.6))
        ax.plot(vs, V_A, "-", lw=1.6, label="V0*sqrt(1-v)  (Model A, QM)")
        ax.plot(vs, V_B, "--", lw=1.6, label="2(1-v)/(2-v)  (Model B)")
        for o in [0.1, 0.25, 0.5, 0.75, 0.9]:
            m = ModelA(o, seed=a.seed)
            vt, v1, v2, v3, V, _ = report_point(m, f"A(o={o:.2f})", a.N, a.ncal, a.tol)
            ax.plot(vt, V, "o", color="tab:blue", ms=6)
        for v in [0.1, 0.25, 0.5, 0.75, 0.9]:
            m = ModelB(v, seed=a.seed)
            vt, v1, v2, v3, V, _ = report_point(m, f"B(v={v:.2f})", a.N, a.ncal, a.tol)
            ax.plot(vt, V, "s", color="tab:red", ms=6)
        ax.set_xlabel("Ze velocity v")
        ax.set_ylabel("measured fringe visibility V")
        ax.set_title("Visibility law: Model A (QM) vs Model B (materialization)")
        ax.legend()
        fig.tight_layout()
        fig.savefig(a.plot, dpi=150)
        print(f"[plot saved: {a.plot}]")


if __name__ == "__main__":
    main()
