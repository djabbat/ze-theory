#!/usr/bin/env python3
"""
ПРЯМОЕ ДОКАЗАТЕЛЬСТВО: H_Ze на пирохлорной решётке → U(1) спиновая жидкость
======================================================================

НЕ через связь с XXZ-моделью.
ПРЯМОЕ измерение g(Γ) из H_Ze через гексагонные корреляторы.

Метод: QMC + измерение ⟨B(0)B(r)⟩ → извлечение g(Γ).

Автор: Jaba Tqemaladze, MD | 2026-07-04
"""
import numpy as np, math, json, sys, os, time
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'pyrochlore'))
try:
    from pyro_lattice import build_pyrochlore
    from ze_qmc_pyro import setup_trotter, wolff_cluster_pyro
    HAS_PYRO = True
except ImportError:
    HAS_PYRO = False
    print("⚠️  pyro_lattice not found")

from numba import njit

LN2 = math.log(2)
V_STAR = 1.0 - LN2
PT_ERROR = (2 - LN2) / 2
ALPHA_EXP = 1 / 137.035999084

@njit
def measure_hex_operator(z, hexagons, M):
    """Прямое измерение ⟨B_p⟩ для всех гексагонов."""
    nh = len(hexagons)
    B_vals = np.zeros(nh)
    for h in range(nh):
        Bh = 1.0
        for tau in range(M):
            prod = 1.0
            for node in hexagons[h]:
                prod *= z[node, tau]
            Bh *= prod ** (1.0 / M)  # геометрическое среднее по tau
        B_vals[h] = Bh
    return np.mean(B_vals), np.std(B_vals) / np.sqrt(nh)

@njit
def measure_hex_correlator_direct(z, hexagons, centers, M, max_dist=12):
    """Измерение ⟨B(0)B(r)⟩ — гексагонный коррелятор."""
    nh = len(hexagons)
    corr = np.zeros(max_dist + 1)
    cnt = np.zeros(max_dist + 1)
    
    for tau in range(M):
        B = np.ones(nh)
        for h in range(nh):
            for node in hexagons[h]:
                B[h] *= z[node, tau]
        
        for i in range(nh):
            for j in range(i+1, nh):
                dx = centers[i,0] - centers[j,0]
                dy = centers[i,1] - centers[j,1]
                dz = centers[i,2] - centers[j,2]
                d = int(math.sqrt(dx*dx+dy*dy+dz*dz) + 0.5)
                if 0 <= d <= max_dist:
                    corr[d] += B[i] * B[j]
                    cnt[d] += 1
    
    for d in range(max_dist + 1):
        if cnt[d] > 0:
            corr[d] /= cnt[d]
    return corr, cnt

def extract_g_from_hex_correlator(corr, distances, L_system):
    """
    Извлечение g из гексагонного коррелятора.
    
    В U(1)-фазе: ⟨B(0)B(r)⟩ ∼ 1/r⁴ на больших расстояниях.
    При r ∼ a: ⟨B(0)B(a)⟩ ≈ 1 − c·g.
    
    Альтернативно: g = (1 − ⟨B(0)B(a)⟩)/c, где c ∼ O(1).
    """
    # Если ice-rule выполняется: corr[d] ≈ 1 для всех d
    if len(corr) >= 3 and np.all(np.abs(corr[1:4] - 1.0) < 0.02):
        return 0.0, 'U(1) phase (plateau, g→0)'
    
    # Экспоненциальный спад
    if corr[1] > 0.02:
        log_c = np.log(np.abs(corr[1:]) + 1e-15)
        d_vals = distances[1:]
        mask = log_c < -0.01
        if mask.any():
            coeffs = np.polyfit(d_vals[mask], log_c[mask], 1)
            xi = -1.0/coeffs[0] if coeffs[0] < 0 else L_system
            # g ∼ 1/ξ для U(1)-жидкости
            g = 1.0 / xi if xi > 0 else 0.0
            return g, f'ξ={xi:.1f}a, g=1/ξ={g:.4f}'
    
    return 0.0, 'no decay detected'

def direct_g_measurement(L=4, Gamma=0.94, M=64, n_thermal=3000, n_meas=500):
    """
    ПРЯМОЕ измерение g из H_Ze на пирохлорной решётке.
    БЕЗ использования данных Hermele et al.
    """
    if not HAS_PYRO:
        return None
    
    print(f"\n{'='*65}")
    print(f"ПРЯМОЕ ИЗМЕРЕНИЕ g(Γ) ИЗ H_Ze: L={L}, Γ={Gamma}")
    print(f"{'='*65}")
    
    t0 = time.time()
    neigh, deg, pos = build_pyrochlore(L)
    N = len(deg)
    
    # Поиск гексагонов
    from numba import njit
    
    @njit
    def find_hex(neigh, deg):
        N = len(deg)
        hexs = []
        for start in range(min(N, 100)):
            paths = [[start]]
            for depth in range(1, 7):
                new = []
                for p in paths:
                    last = p[-1]
                    for k in range(deg[last]):
                        nb = neigh[last, k]
                        if depth < 6:
                            if nb not in p: new.append(p + [nb])
                        else:
                            if nb == start:
                                dup = False
                                for h in hexs:
                                    if len(h)==6:
                                        for s in range(6):
                                            if all(p[i]==h[(i+s)%6] for i in range(6)):
                                                dup = True; break
                                            if all(p[i]==h[(s-i)%6] for i in range(6)):
                                                dup = True; break
                                if not dup: hexs.append(p.copy())
                paths = new
                if not paths: break
        return hexs
    
    hexagons_raw = find_hex(neigh, deg)
    hexagons = np.zeros((len(hexagons_raw), 6), dtype=np.int32)
    for hi, h in enumerate(hexagons_raw):
        for k, n in enumerate(h): hexagons[hi, k] = n
    
    centers = np.zeros((len(hexagons), 3))
    for h in range(len(hexagons)):
        cx = cy = cz = 0.0
        for node in hexagons[h]: cx += pos[node,0]; cy += pos[node,1]; cz += pos[node,2]
        centers[h] = [cx/6, cy/6, cz/6]
    
    print(f"  N={N}, гексагонов={len(hexagons)}")
    
    # QMC
    Ks, Kt = setup_trotter(1.0, Gamma, M)
    np.random.seed(42)
    z = np.random.choice(np.array([-1.0, 1.0]), size=(N, M))
    
    print(f"  Термализация ({n_thermal} шагов)...")
    for step in range(n_thermal):
        wolff_cluster_pyro(z, neigh, deg, M, Ks, Kt)
    
    print(f"  Измерения ({n_meas} × 5 шагов)...")
    B_means = np.zeros(n_meas)
    max_dist = min(12, L*3)
    all_corr = np.zeros((n_meas, max_dist+1))
    
    for k in range(n_meas):
        for _ in range(5): wolff_cluster_pyro(z, neigh, deg, M, Ks, Kt)
        B_means[k], _ = measure_hex_operator(z, hexagons, M)
        corr, _ = measure_hex_correlator_direct(z, hexagons, centers, M, max_dist)
        all_corr[k] = corr[:max_dist+1]
    
    mean_B = np.mean(B_means)
    mean_corr = np.mean(all_corr, axis=0)
    
    # Извлечение g
    distances = np.arange(max_dist+1)
    g_direct, desc = extract_g_from_hex_correlator(mean_corr, distances, L)
    
    elapsed = time.time() - t0
    alpha_direct = PT_ERROR * abs(g_direct) / (4 * math.pi) if abs(g_direct) > 1e-12 else 0
    
    print(f"\n  [{elapsed:.0f}s] Результаты:")
    print(f"  ⟨B⟩ = {mean_B:.6f}")
    print(f"  ⟨B(0)B(1)⟩ = {mean_corr[1]:.6f}")
    print(f"  g(Γ={Gamma}) = {g_direct:.6f} ({desc})")
    print(f"  α(g) = {alpha_direct:.8f}" + (f" (1/{1/alpha_direct:.1f})" if alpha_direct > 0 else ""))
    
    return {
        'L': L, 'Gamma': Gamma, 'M': M,
        'mean_B': float(mean_B),
        'corr_1': float(mean_corr[1]),
        'g_direct': float(g_direct),
        'alpha_direct': float(alpha_direct),
        'alpha_exp': ALPHA_EXP,
        'g_desc': desc, 'time': elapsed,
        'method': 'DIRECT measurement from H_Ze, NO Hermele data used'
    }

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description='Direct g(Γ) measurement from H_Ze')
    p.add_argument('--L', type=int, default=4)
    p.add_argument('--Gamma', type=float, default=0.94)
    p.add_argument('--save', type=str, default=None)
    args = p.parse_args()
    
    r = direct_g_measurement(L=args.L, Gamma=args.Gamma)
    if r and args.save:
        with open(args.save, 'w') as f: json.dump(r, f, indent=2)
        print(f"Saved: {args.save}")
    elif r:
        print(f"\ng={r['g_direct']:.6f}, α={r['alpha_direct']:.8f}")
