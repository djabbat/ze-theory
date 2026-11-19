#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Decoherence = phase desynchronization = growing prediction error (Ze).

Formalization (FREQUENCY_ERROR_Ze.md, checklist item 2):
  reality  p(x,t)  = env(x) (1 + V(t) cos(delta(x))) / Z,
                    V(t) = V0 exp(-delta_phi(t)^2 / 2)     (Gaussian dephasing)
  model    p_mod(x)= env(x) (1 + cos(delta(x))) / Z        (assumes full coherence)
  error    I_n(t)  = D_KL(p_real(t) || p_mod)

  delta_phi(t) = delta_E * t / hbar   (phase desynchronization from energy spread)
  tau_dec = hbar / delta_E            (timescale on which prediction error grows)

Ze reading: a Ze-system whose model keeps assuming full coherence while the
phases desynchronize accumulates prediction error I_n(t); decoherence IS this
growing error, and its timescale is set by the frequency spread (energy-time
uncertainty). This connects to the CHSH-entropy law (injected entropy H ~
phase/error noise degrades correlations).

Run: python3 decoherence_error.py
Only numpy is required.
"""

import math
import numpy as np

LAMBDA_L = 1.0
SLIT_WIDTH = 0.2
SLIT_SEP = 0.5
XMAX = 4.9
NPIX = 401


def grid():
    return np.linspace(-XMAX, XMAX, NPIX)


def envelope(x):
    u = math.pi * SLIT_WIDTH * x / LAMBDA_L
    out = np.ones_like(x)
    m = np.abs(u) > 1e-12
    out[m] = (np.sin(u[m]) / u[m]) ** 2
    return out


def delta_phase(x):
    return 2.0 * math.pi * SLIT_SEP * x / LAMBDA_L


def pdf(x, V):
    p = envelope(x) * (1.0 + V * np.cos(delta_phase(x)))
    return p / p.sum()


def kl(p, q):
    m = p > 1e-15
    return float(np.sum(p[m] * np.log(p[m] / q[m])))


def test_decoherence_error():
    print("[decoherence = phase desync = growing prediction error]")
    x = grid()
    p_mod = pdf(x, 1.0)          # model: assumes full coherence (V=1)
    # energy spread in natural units hbar = 1
    results = []
    for dE in (0.1, 0.2, 0.5):
        tau = 1.0 / dE           # hbar/dE
        ts = np.linspace(0.0, 5.0 * tau, 60)
        In = np.array([kl(pdf(x, math.exp(-(dE * t) ** 2 / 2.0)), p_mod)
                       for t in ts])
        # small-t growth ~ c t^2 (Gaussian decoherence signature)
        small = ts < 0.3 * tau
        c = float(np.polyfit(ts[small], In[small], 2)[0])
        I_inf = float(In[-1])
        idx = int(np.where(In > 0.5 * I_inf)[0][0])
        t_half = float(ts[idx])
        results.append((dE, tau, I_inf, c, t_half / tau))
        print(f"   dE={dE:.1f} tau_dec=hbar/dE={tau:.2f}: "
              f"I_n max={I_inf:.3f}, I_n~c*t^2 c={c:.3f}, "
              f"t(I_n=half)/tau={t_half / tau:.3f}")
    print("   => I_n(t) grows from 0 (model perfect at t=0) to a plateau;")
    print("      small-t growth is quadratic (Gaussian dephasing); the")
    print("      timescale is tau_dec = hbar/dE (energy-time uncertainty).")
    print("   Ze: decoherence = accumulated prediction error of a model that")
    print("      ignores phase desynchronization.")
    # PASS criteria: (a) error grows to a plateau; (b) t_half/tau scale-invariant
    # (same decoherence timescale for all dE); (c) c ~ dE^2 (Gaussian dephasing).
    ratios = [r[4] for r in results]
    scale_inv = (max(ratios) - min(ratios)) / max(ratios) < 0.05
    c_over_dE2 = [r[3] / r[0] ** 2 for r in results]
    quad = (max(c_over_dE2) - min(c_over_dE2)) / max(c_over_dE2) < 0.15
    ok = all(r[2] > 0 for r in results) and scale_inv and quad
    print(f"   scale-invariance t_half/tau: {'OK' if scale_inv else 'FAIL'}; "
          f"c~dE^2: {'OK' if quad else 'FAIL'}")
    print(f"   VERDICT: {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    test_decoherence_error()


if __name__ == "__main__":
    main()
