"""
Анализ результатов Ze-Монте-Карло: поиск v* и фазовая диаграмма.
Запускать после ze_mc.py с генерацией results_*.json
"""
import json
import numpy as np
from pathlib import Path


def load_results(filepath):
    with open(filepath) as f:
        return json.load(f)


def find_v_star_crossing(results_list):
    """Найти параметры, где |v| ≈ v* = 1 − ln 2."""
    v_star = 1.0 - np.log(2.0)
    print(f"v* = {v_star:.6f}")
    print()
    
    for r in results_list:
        v = r.get("magnetization_abs", abs(r.get("mean_magnetization", 0)))
        diff = abs(v - v_star)
        p = r["params"]
        print(f"T={p['T']:.3f} J_s={p['J_s']:.3f} h={p['h']:.3f}  "
              f"|v|={v:.4f}  ξ_t={r.get('xi_t', 0):.2f}  "
              f"Δv={diff:.4f} {'← БЛИЗКО!' if diff < 0.05 else ''}")


# ============================================================
# Интерактивный анализ
# ============================================================
if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        files = sys.argv[1:]
    else:
        files = list(Path(".").glob("results_*.json"))
    
    all_results = [load_r(f) for f in files]
    find_v_star_crossing(all_results)
