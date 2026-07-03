"""
3+1d Монте-Карло Ze-модели с измерением петель Вильсона.

H = -J_t Σ z(x,y,z,t)z(x,y,z,t+1) - J_s Σ [z_x+...+z_y+...+z_z+...] - h Σ z

Петли Вильсона: W(C) = ⟨Π_{(i,j)∈C} z_i z_j⟩
- Perimeter law → деконфайнмент (U(1)-фаза)
- Area law → конфайнмент
"""

import numpy as np
from numba import njit
from datetime import datetime
import json

# ============================================================
@njit
def energy_4d(z, L, Lt, J_t, J_s, h):
    """Полная энергия конфигурации."""
    E = 0.0
    for x in range(L):
        for y in range(L):
            for zz in range(L):
                for t in range(Lt):
                    # поле
                    E -= h * z[x, y, zz, t]
                    # временная связь
                    E -= J_t * z[x, y, zz, t] * z[x, y, zz, (t+1) % Lt]
                    # пространственные связи
                    E -= J_s * z[x, y, zz, t] * z[(x+1) % L, y, zz, t]
                    E -= J_s * z[x, y, zz, t] * z[x, (y+1) % L, zz, t]
                    E -= J_s * z[x, y, zz, t] * z[x, y, (zz+1) % L, t]
    return E


@njit
def metropolis_sweep_4d(z, L, Lt, J_t, J_s, h, beta, n_sweeps=1):
    """Свипы Метрополиса."""
    total = L * L * L * Lt
    acc = 0
    for _ in range(n_sweeps):
        for _ in range(total):
            idx = np.random.randint(0, total)
            t = idx % Lt; idx //= Lt
            zz = idx % L; idx //= L
            y = idx % L; x = idx // L
            
            # локальная энергия
            e_old = -h * z[x, y, zz, t]
            e_old -= 2*J_t * z[x, y, zz, t] * z[x, y, zz, (t+1) % Lt]
            e_old -= 2*J_s * z[x, y, zz, t] * z[(x+1) % L, y, zz, t]
            e_old -= 2*J_s * z[x, y, zz, t] * z[x, (y+1) % L, zz, t]
            e_old -= 2*J_s * z[x, y, zz, t] * z[x, y, (zz+1) % L, t]
            
            z[x, y, zz, t] = -z[x, y, zz, t]
            
            e_new = -h * z[x, y, zz, t]
            e_new -= 2*J_t * z[x, y, zz, t] * z[x, y, zz, (t+1) % Lt]
            e_new -= 2*J_s * z[x, y, zz, t] * z[(x+1) % L, y, zz, t]
            e_new -= 2*J_s * z[x, y, zz, t] * z[x, (y+1) % L, zz, t]
            e_new -= 2*J_s * z[x, y, zz, t] * z[x, y, (zz+1) % L, t]
            
            if e_new <= e_old or np.random.random() < np.exp(-beta * (e_new - e_old)):
                acc += 1
            else:
                z[x, y, zz, t] = -z[x, y, zz, t]
    return acc / total


@njit
def measure_wilson_loops(z, L, Lt, max_size=3):
    """Измерение петель Вильсона. Возвращает массивы."""
    n_loops = max_size * max_size
    areas = np.zeros(n_loops)
    perimeters = np.zeros(n_loops)
    W_values = np.zeros(n_loops)
    idx = 0
    
    for R in range(1, max_size + 1):
        for T_val in range(1, max_size + 1):
            w_sum = 0.0
            count = 0
            
            for x in range(L - R):
                for y in range(L):
                    for zz in range(L):
                        for t in range(Lt - T_val):
                            # x-t loop
                            loop = 1.0
                            for dx in range(R):
                                loop *= z[x+dx, y, zz, t]
                            for dt in range(T_val):
                                loop *= z[x+R, y, zz, t+dt]
                            for dx in range(R):
                                loop *= z[x+R-dx, y, zz, t+T_val]
                            for dt in range(T_val):
                                loop *= z[x, y, zz, t+T_val-dt]
                            w_sum += loop
                            count += 1
            
            if count > 0:
                areas[idx] = R * T_val
                perimeters[idx] = 2 * (R + T_val)
                W_values[idx] = w_sum / count
                idx += 1
    
    return areas[:idx], perimeters[:idx], W_values[:idx]


def run_4d_simulation(L=4, Lt=8, J_t=1.0, J_s=0.3, h=0.0, T=2.0,
                     n_thermal=2000, n_samples=5000, sample_int=10):
    """Запуск 4D симуляции."""
    beta = 1.0 / T
    np.random.seed(42)
    z = np.random.choice([-1, 1], size=(L, L, L, Lt)).astype(np.float64)
    
    print(f"3+1d MC: L={L}, Lt={Lt}, J_t={J_t}, J_s={J_s}, T={T}")
    N = L*L*L*Lt
    print(f"  Всего узлов: {N}")
    
    # термализация
    print(f"  Термализация ({n_thermal} шагов)...")
    metropolis_sweep_4d(z, L, Lt, J_t, J_s, h, beta, n_sweeps=n_thermal)
    
    # измерения
    print(f"  Измерения ({n_samples} шагов)...")
    n_meas = n_samples // sample_int
    energies = np.zeros(n_meas)
    magnetizations = np.zeros(n_meas)
    
    for i in range(n_samples):
        metropolis_sweep_4d(z, L, Lt, J_t, J_s, h, beta, n_sweeps=1)
        if i % sample_int == 0:
            idx = i // sample_int
            energies[idx] = energy_4d(z, L, Lt, J_t, J_s, h) / N
            magnetizations[idx] = abs(np.mean(z))
    
    # петли Вильсона
    print(f"  Измерение петель Вильсона...")
    areas, perimeters, W_vals = measure_wilson_loops(z, L, Lt, max_size=min(3, L//2))
    
    # анализ шлейфа Вильсона
    v_star = 1.0 - np.log(2.0)
    results = {
        "params": {"L": L, "Lt": Lt, "J_t": J_t, "J_s": J_s, "h": h, "T": T},
        "E_per_node": float(np.mean(energies)),
        "E_std": float(np.std(energies)),
        "v_abs": float(np.mean(magnetizations)),
        "v_std": float(np.std(magnetizations)),
        "v_star": v_star,
        "delta_v": float(abs(np.mean(magnetizations) - v_star)),
    }
    
    # анализ фазового состояния
    if len(W_vals) > 0:
        mask = (np.abs(W_vals) > 1e-10) & (areas > 0)
        if mask.sum() >= 2:
            coeffs = np.polyfit(areas[mask], np.log(np.abs(W_vals[mask])), 1)
            results["string_tension"] = float(-coeffs[0])
            results["phase"] = "confinement" if results["string_tension"] > 0.01 else "possible deconfinement"
        # также фит по периметру для сравнения
        if mask.sum() >= 2:
            coeffs_p = np.polyfit(perimeters[mask], np.log(np.abs(W_vals[mask])), 1)
            results["perimeter_coeff"] = float(-coeffs_p[0])
    
    print(f"\n  ⟨E⟩/N = {results['E_per_node']:.4f}")
    print(f"  |v| = {results['v_abs']:.4f} (v*={v_star:.4f})")
    if "string_tension" in results:
        print(f"  string tension σ = {results['string_tension']:.4f}")
        print(f"  фаза: {results['phase']}")
    
    return results


if __name__ == "__main__":
    # сканирование по J_s
    print("Поиск U(1)-фазы в 3+1d Ze-модели\n")
    
    for J_s in [0.1, 0.2, 0.3, 0.5, 0.7]:
        for T in [2.0, 2.5, 3.0]:
            r = run_4d_simulation(L=4, Lt=6, J_t=1.0, J_s=J_s, T=T,
                                 n_thermal=1000, n_samples=3000)
            phase = r.get("phase", "?")
            sigma = r.get("string_tension", float('nan'))
            print(f"  J_s={J_s} T={T}: v={r['v_abs']:.3f} σ={sigma:.4f} {phase}")
            print()
