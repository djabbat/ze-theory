#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compton frequency from the Ze equation of motion — numerical check.

Formalization (ENERGY_MASS_FREQUENCY_Ze.md, checklist item 1):
  * The Ze EOM (Schrodinger, derived from entropic dynamics, manuscript §4.5)
    admits stationary states  psi(x,t) = sqrt(rho)(x) e^{-i E t / hbar}.
  * The phase rotates at  omega = E / hbar = 2 pi nu.
  * For a particle at rest  E = m c^2  (mass-energy identification)
      => nu_c = m c^2 / h   (Compton frequency).
  * In Ze: the rest clock ticks at nu_c; it runs error-free (stationary,
    N_S = 0, grad I_n = 0); proper time advances at the phase rate
    tau = N_T (all events are T/stasis). The frequency-energy relation
    nu = E/h is DERIVED from the EOM; m and E = m c^2 are inputs.

Check: evolve a bound stationary state with the unitary Crank-Nicolson
integrator (the Ze EOM) and measure the phase rotation rate of psi at a
fixed point; verify it equals E / hbar = 2 pi nu.

For the harmonic oscillator (omega_osc), the ground state has
E = hbar omega_osc / 2, so the phase rotates at omega = omega_osc / 2
and the frequency is nu = omega_osc / (4 pi).

Run: python3 compton_frequency.py
Only numpy is required.
"""

import math
import numpy as np

# ----------------------------------------------------------------------
# Crank-Nicolson unitary integrator (Ze EOM), same as ed_two_slit.py
# ----------------------------------------------------------------------
def laplacian(psi, h):
    out = np.zeros_like(psi)
    out[1:-1] = (psi[:-2] - 2 * psi[1:-1] + psi[2:]) / h ** 2
    return out


def thomas(a, b, c, d):
    n = len(d)
    cp = np.zeros(n - 1, dtype=complex)
    dp = d.copy().astype(complex)
    cp[0] = c[0] / b[0]
    dp[0] = dp[0] / b[0]
    for i in range(1, n):
        m = b[i] - a[i - 1] * cp[i - 1]
        if i < n - 1:
            cp[i] = c[i] / m
        dp[i] = (d[i] - a[i - 1] * dp[i - 1]) / m
    x = np.zeros(n, dtype=complex)
    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


def cn_step(psi, V, dt, h):
    n = len(psi)
    main = np.ones(n) / h ** 2 + V
    off = -0.5 * np.ones(n - 1) / h ** 2
    rhs = psi - 0.5j * dt * (-0.5 * laplacian(psi, h) + V * psi)
    a = 0.5j * dt * off
    b = 1.0 + 0.5j * dt * main
    c = 0.5j * dt * off
    return thomas(a, b, c, rhs)


def energy(psi, V, h):
    ke = -0.5 * np.real(np.vdot(psi, laplacian(psi, h)))
    pe = np.real(np.vdot(psi, V * psi))
    return ke + pe


def main():
    # harmonic oscillator: ground state E = hbar*omega/2 (hbar=1)
    omega = 1.0
    x = np.linspace(-12.0, 12.0, 801)
    h = x[1] - x[0]
    V = 0.5 * omega ** 2 * x ** 2
    psi = np.exp(-omega * x ** 2 / 2.0)
    psi /= np.sqrt(np.vdot(psi, psi))
    E = energy(psi, V, h)
    E_exact = 0.5 * omega          # ground-state energy (hbar=1)
    # evolve and measure the phase of psi at x=0 over time
    dt = 0.005
    n_steps = 1200
    phases = []
    ts = []
    for k in range(n_steps + 1):
        phases.append(math.atan2(psi[psi.size // 2].imag, psi[psi.size // 2].real))
        ts.append(k * dt)
        psi = cn_step(psi, V, dt, h)
    phases = np.unwrap(np.asarray(phases))
    ts = np.asarray(ts)
    # linear fit of phase vs time -> omega_meas
    omega_meas = float(np.polyfit(ts, phases, 1)[0])
    print("[Compton frequency from the Ze EOM]")
    print(f"   harmonic oscillator: E_num={E:.6f} (exact {E_exact:.6f})")
    print(f"   measured phase rate omega_meas = {omega_meas:.6f} "
          f"(phase of e^(-iEt/hbar) decreases -> negative rate)")
    print(f"   |omega_meas| = E/hbar = {abs(omega_meas):.6f} (expected {E_exact:.6f})")
    print(f"   nu = |omega|/(2 pi) = {abs(omega_meas) / (2 * math.pi):.6f}")
    ok = abs(abs(omega_meas) - E_exact) / E_exact < 1e-3
    print(f"   => phase rotation rate = E/hbar  {'PASS' if ok else 'FAIL'}")
    print("   => for a particle at rest (E = m c^2): nu_c = m c^2 / h")
    print("      (electron: nu_c ~ 1.24e20 Hz) -- the rest clock rate from")
    print("      the Ze equation of motion.")
    print(f"   VERDICT: {'PASS' if ok else 'FAIL'}")
    return ok


if __name__ == "__main__":
    main()
