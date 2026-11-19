#!/usr/bin/env python3
"""S1 seam: quantized phase. Sensitivity of a fringe measurement to the Ze register step.

Model. The phase difference between the two alternatives is a count difference:
    delta(f) = 2*pi*f,   f = fringe order (accumulated phase in units of 2*pi)
If the register advances in discrete steps and N is the number of steps per
full 2*pi of phase, the realizable phase is the staircase
    f_N = round(N*f)/N,      p_N(f) = cos^2(pi*f_N),
against the continuous p(f) = cos^2(pi*f).

Outputs: RMS/max deviation, visibility of the staircase, and the step count
N required before the deviation drops below a stated accuracy. The step length
is then l = lambda/N, which is the experimental bound S1 sets.
"""
import numpy as np

ACC = (0.1, 0.01, 0.001, 1e-4)
Ns = [1, 2, 3, 5, 10, 20, 50, 100, 300, 1000, 3000, 1e4, 1e5, 1e6]

f = np.linspace(0.0, 1.0, 200001)
p = np.cos(np.pi * f) ** 2

print("N steps/2pi |  RMS dev  |  max dev  | visibility | fringe-pos err (frac of order)")
print("-" * 88)
rows = []
for N in Ns:
    fq = np.round(N * f) / N
    pq = np.cos(np.pi * fq) ** 2
    rms = float(np.sqrt(np.mean((pq - p) ** 2)))
    mx = float(np.max(np.abs(pq - p)))
    vis = float(pq.max() - pq.min())
    rows.append((N, rms, mx, vis))

for N, rms, mx, vis in rows:
    print(f"{int(N):>10} | {rms:9.5f} | {mx:9.5f} | {vis:10.5f} | -")

print()
print("Step count required for a given accuracy of the fringe shape (max deviation):")
for a in ACC:
    Nreq = next((N for N, rms, mx, vis in rows if mx <= a), None)
    # analytic expectation: max deviation ~ pi/N  =>  N ~ pi/a
    print(f"  max dev <= {a:g}:  N >= {int(Nreq) if Nreq else '>1e6'}   (analytic ~ pi/a = {np.pi/a:.3g})")

print()
lam = {"100 keV electron": 3.70e-12, "thermal C60 (900 m/s)": 4.6e-12,
       "cold atom (Rb, 1 mm/s)": 6.4e-7}
print("Bound on the register step  l = lambda / N  (for max dev <= 1%):")
N1 = next((N for N, rms, mx, vis in rows if mx <= 0.01), 3.2e2)
for name, L in lam.items():
    print(f"  {name:26s} lambda = {L:9.3e} m ->  l <= {L/N1:9.3e} m")
print()
print(f"Planck length 1.6e-35 m would give N ~ lambda/l_P ~ {3.7e-12/1.6e-35:.1e}"
      f" steps per wavelength: the seam is invisible at Planck-scale steps.")
print("Inversion: any measured max deviation `d` at precision `d` sets N >= pi/d, hence l <= lambda*d/pi.")
