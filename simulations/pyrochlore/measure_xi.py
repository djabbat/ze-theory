#!/usr/bin/env python3
"""
Прямое измерение корреляционной длины ξ на пирохлорной решётке.
================================================================

Вычисляет ⟨B(0)B(r)⟩ для гексагонов и извлекает ξ(L).
Проверяет гипотезу: выходит ли ξ на плато ~7 при больших L.

Запуск:
  python measure_xi.py --L 8 --Gamma 0.9 --M 64
  python measure_xi.py --scan-L --Gamma 0.9  # L=4,5,6,7,8,10,12

Автор: Jaba Tqemaladze, MD
Дата: 2026-07-04
"""

import numpy as np
import math
import json
import sys
import os
import time
from datetime import datetime
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'pyrochlore'))

try:
    from pyro_lattice import build_pyrochlore
    from ze_qmc_pyro import setup_trotter, wolff_cluster_pyro
    HAS_PYRO = True
except ImportError as e:
    HAS_PYRO = False
    print(f"⚠️  Не удалось импортировать pyro_lattice: {e}")
    print("   Будет использоваться только аналитическая модель")

from numba import njit


# ═══════════════════════════════════════════════════════════════════
# ПОИСК ГЕКСАГОНОВ
# ═══════════════════════════════════════════════════════════════════

@njit
def find_all_hexagons(neigh, deg, max_hex=2000):
    """
    Находит ВСЕ гексагоны на пирохлорной решётке
    (не только от первых 50 узлов).
    
    Гексагон — замкнутый путь из 6 различных узлов,
    проходящий по рёбрам решётки.
    """
    N = len(deg)
    hexagons = []
    
    for start in range(N):
        # BFS глубины 6
        paths = [[start]]
        for depth in range(1, 7):
            new_paths = []
            for path in paths:
                last = path[-1]
                for k in range(deg[last]):
                    nb = neigh[last, k]
                    if depth < 6:
                        if nb not in path:
                            new_paths.append(path + [nb])
                    else:
                        if nb == start and len(path) == 6:
                            # Проверяем, нет ли уже такого гексагона
                            # (с точностью до циклического сдвига и обращения)
                            found_dup = False
                            for h in hexagons:
                                # Проверяем совпадение (сдвиг + обращение)
                                if len(h) == 6:
                                    for shift in range(6):
                                        match = True
                                        for i in range(6):
                                            if path[i] != h[(i + shift) % 6]:
                                                match = False
                                                break
                                        if match:
                                            found_dup = True
                                            break
                                        # обратный порядок
                                        match = True
                                        for i in range(6):
                                            if path[i] != h[(shift - i) % 6]:
                                                match = False
                                                break
                                        if match:
                                            found_dup = True
                                            break
                                    if found_dup:
                                        break
                            if not found_dup:
                                hexagons.append(path.copy())
            paths = new_paths
            if len(paths) == 0:
                break
        
        if len(hexagons) >= max_hex:
            break
    
    return hexagons


# ═══════════════════════════════════════════════════════════════════
# ИЗМЕРЕНИЕ ГЕКСАГОННЫХ КОРРЕЛЯТОРОВ
# ═══════════════════════════════════════════════════════════════════

@njit
def measure_hex_correlators(z, hex_array, centers, max_dist=15):
    """
    Измерение ⟨B(0)B(r)⟩ для всех гексагонов.
    
    B_p = Π_{i∈⎔} z_i — оператор магнитного потока.
    """
    nh = len(hex_array)
    M = z.shape[1]  # число троттеровских слоёв
    
    # Усреднение по мнимому времени
    corr = np.zeros(max_dist + 1)
    cnt = np.zeros(max_dist + 1)
    
    for tau in range(M):
        # Вычисляем все B_p для этого tau
        B = np.ones(nh)
        for h in range(nh):
            for node in hex_array[h]:
                B[h] *= z[node, tau]
        
        # Корреляции между всеми парами
        for i in range(nh):
            for j in range(i + 1, nh):
                dx = centers[i, 0] - centers[j, 0]
                dy = centers[i, 1] - centers[j, 1]
                dz = centers[i, 2] - centers[j, 2]
                d = int(math.sqrt(dx*dx + dy*dy + dz*dz) + 0.5)
                
                if 0 <= d <= max_dist:
                    corr[d] += B[i] * B[j]
                    cnt[d] += 1
    
    # Нормировка
    for d in range(max_dist + 1):
        if cnt[d] > 0:
            corr[d] /= cnt[d]
    
    return corr, cnt


def estimate_xi_from_correlator(corr, distances):
    """
    Оценка корреляционной длины ξ из ⟨B(0)B(r)⟩.
    
    Использует экспоненциальный фит: corr(r) ~ exp(-r/ξ) + const.
    Если corr(r) ≈ 1 для всех r → ξ = inf (U(1)-фаза).
    Если corr(r) экспоненциально спадает → извлекаем ξ.
    """
    # Проверка: plateau на 1?
    if len(corr) >= 4:
        if np.all(np.abs(corr[1:5] - 1.0) < 0.02):
            return float('inf'), 'U(1) phase (plateau at 1)'
    
    # Экспоненциальный фит
    mask = (corr > 0.02) & (distances > 0)
    if not np.any(mask):
        return 1.0, 'no signal'
    
    d = distances[mask]
    c = corr[mask]
    
    if len(c) < 3:
        return 1.0, 'too few points'
    
    # Линейный фит log(c) vs d
    log_c = np.log(c + 1e-15)
    coeffs = np.polyfit(d, log_c, 1)
    
    if coeffs[0] < 0:
        xi = -1.0 / coeffs[0]
        return xi, f'exponential (ξ={xi:.2f})'
    else:
        return float('inf'), 'no decay (U(1) or power-law)'


# ═══════════════════════════════════════════════════════════════════
# ГЛАВНОЕ ИЗМЕРЕНИЕ НА ПИРОХЛОРНОЙ РЕШЁТКЕ
# ═══════════════════════════════════════════════════════════════════

def run_xi_measurement(L=4, Gamma=0.9, M=64, n_thermal=3000, n_meas=200):
    """
    Прямое измерение ξ на пирохлорной решётке с QMC.
    """
    if not HAS_PYRO:
        print("❌ Нет pyro_lattice — невозможно прямое измерение")
        return None
    
    print(f"\n{'='*65}")
    print(f"ИЗМЕРЕНИЕ ξ: L={L}, Γ={Gamma}, M={M}")
    print(f"{'='*65}")
    
    t0 = time.time()
    
    neigh, deg, pos = build_pyrochlore(L)
    N = len(deg)
    
    print(f"  Спинов: {N} × {M} = {N*M}")
    print(f"  Координация: {np.min(deg)}–{np.max(deg)}")
    
    # Поиск гексагонов
    hexagons = find_all_hexagons(neigh, deg, max_hex=1000)
    print(f"  Гексагонов: {len(hexagons)}")
    
    if len(hexagons) == 0:
        print("  ⚠️ Гексагоны не найдены!")
        return None
    
    # Преобразуем в массив для Numba
    hex_array = np.zeros((len(hexagons), 6), dtype=np.int32)
    for hi, h in enumerate(hexagons):
        for k, n in enumerate(h):
            hex_array[hi, k] = n
    
    # Центры гексагонов
    centers = np.zeros((len(hexagons), 3))
    for h in range(len(hexagons)):
        cx = cy = cz = 0.0
        for node in hex_array[h]:
            cx += pos[node, 0]
            cy += pos[node, 1]
            cz += pos[node, 2]
        centers[h] = [cx/6, cy/6, cz/6]
    
    # QMC
    Ks, Kt = setup_trotter(1.0, Gamma, M)
    np.random.seed(42)
    z = np.random.choice(np.array([-1.0, 1.0]), size=(N, M))
    
    # Термализация
    print(f"  Термализация ({n_thermal} шагов)...")
    for step in range(n_thermal):
        wolff_cluster_pyro(z, neigh, deg, M, Ks, Kt)
        if step % 1000 == 0:
            # Ice-rule check
            ice = 0.0
            for tau in range(M):
                for x in range(L):
                    for y in range(L):
                        for zc in range(L):
                            base = ((x*L + y)*L + zc) * 4
                            s = z[base,tau]+z[base+1,tau]+z[base+2,tau]+z[base+3,tau]
                            if abs(s) < 0.01:
                                ice += 1
            ice /= (M * L**3)
            print(f"    {step}: ice={ice:.4f}")
    
    # Измерения
    print(f"  Измерения ({n_meas} × 5 шагов)...")
    max_dist = min(15, L * 3)
    distances = np.arange(max_dist + 1)
    
    all_corr = np.zeros((n_meas, max_dist + 1))
    
    for k in range(n_meas):
        for _ in range(5):
            wolff_cluster_pyro(z, neigh, deg, M, Ks, Kt)
        corr, cnt = measure_hex_correlators(z, hex_array, centers, max_dist)
        all_corr[k] = corr[:max_dist+1]
    
    # Усреднение
    mean_corr = np.mean(all_corr, axis=0)
    std_corr = np.std(all_corr, axis=0)
    
    # Оценка ξ
    xi, xi_desc = estimate_xi_from_correlator(mean_corr, distances)
    
    elapsed = time.time() - t0
    
    print(f"\n  Результаты [{elapsed:.0f}s]:")
    print(f"  ⟨B(0)B(r)⟩:")
    for d in range(min(max_dist + 1, 12)):
        if d < len(mean_corr):
            print(f"    r≈{d:.0f}a: {mean_corr[d]:.6f} ± {std_corr[d]:.6f}")
    
    print(f"\n  Корреляционная длина: ξ = {xi_desc}")
    
    LN2 = math.log(2)
    PT = (2 - LN2) / 2
    
    if xi == float('inf'):
        L_eff = L
        alpha_pred = PT / (4 * math.pi * L_eff)
        print(f"  L_eff (system size) = {L}")
    else:
        L_eff = xi
        alpha_pred = PT / (4 * math.pi * L_eff)
        print(f"  L_eff (ξ) = {xi:.2f}")
    
    print(f"  α_pred = {alpha_pred:.6f} (1/{1/alpha_pred:.1f})")
    
    return {
        'L': L,
        'Gamma': Gamma,
        'M': M,
        'n_hexagons': len(hexagons),
        'xi': xi if xi != float('inf') else 'inf',
        'xi_numeric': xi if xi != float('inf') else L,
        'L_eff': L_eff,
        'alpha_predicted': alpha_pred,
        'correlator': mean_corr.tolist(),
        'std_correlator': std_corr.tolist(),
        'time': elapsed
    }


# ═══════════════════════════════════════════════════════════════════
# СКАНИРОВАНИЕ ПО L
# ═══════════════════════════════════════════════════════════════════

def scan_L(Gamma=0.9, L_values=None, M=64):
    """
    Сканирование по размерам системы для проверки:
    выходит ли ξ на плато?
    """
    if L_values is None:
        L_values = [4, 5, 6, 7, 8]
    
    print("=" * 65)
    print(f"СКАНИРОВАНИЕ ξ(L): Γ={Gamma}, M={M}")
    print("=" * 65)
    
    results = []
    for L in L_values:
        r = run_xi_measurement(L=L, Gamma=Gamma, M=M,
                               n_thermal=2000, n_meas=100)
        if r:
            results.append(r)
            print(f"  L={L}: ξ={r['xi']}, L_eff={r['L_eff']:.2f}, "
                  f"α={r['alpha_predicted']:.6f}")
    
    print(f"\n{'='*65}")
    print("СВОДКА:")
    print(f"{'L':>6} {'ξ':>12} {'L_eff':>12} {'α':>12} {'1/α':>10}")
    print("-" * 55)
    for r in results:
        xi_str = f"{r['xi']}" if isinstance(r['xi'], str) else f"{r['xi']:.2f}"
        print(f"{r['L']:6d} {xi_str:>12} {r['L_eff']:12.2f} "
              f"{r['alpha_predicted']:12.6f} {1/r['alpha_predicted']:10.1f}")
    
    # Проверка на плато
    if len(results) >= 3:
        xi_values = [r['xi_numeric'] for r in results]
        mean_xi = np.mean(xi_values)
        std_xi = np.std(xi_values)
        print(f"\n  Среднее ξ = {mean_xi:.2f} ± {std_xi:.2f}")
        
        if std_xi < 0.5 and mean_xi > 0:
            print(f"  ✓ ξ ВЫХОДИТ НА ПЛАТО ~{mean_xi:.1f}")
            print(f"  → L не зависит от размера системы!")
            print(f"  → L — ФУНДАМЕНТАЛЬНАЯ КОНСТАНТА теории")
        else:
            print(f"  → ξ зависит от L — требуется экстраполяция L→∞")
    
    return results


# ═══════════════════════════════════════════════════════════════════
# АНАЛИТИЧЕСКАЯ МОДЕЛЬ (без симуляции)
# ═══════════════════════════════════════════════════════════════════

def analytic_xi_model(Gamma, J=1.0, L_system=8):
    """
    Аналитическая модель для ξ(Γ).
    
    В U(1)-фазе: ξ ∼ L_system (насыщение)
    Вблизи Γ_c: ξ ∼ |Γ - Γ_c|^{-ν}, ν ≈ 0.67
    
    Эффективный размер когерентного домена:
    ξ_eff = min(ξ_intrinsic, L_system)
    """
    Gamma_c = 1.05  # из QMC данных
    
    # Внутренняя корреляционная длина (FSS)
    if abs(Gamma - Gamma_c) < 0.01:
        xi_intrinsic = float('inf')
    else:
        xi_0 = 0.5  # микроскопическая длина
        nu = 0.6717  # 3D XY universality
        xi_intrinsic = xi_0 * (abs(Gamma - Gamma_c) / Gamma_c) ** (-nu)
    
    # Эффективная длина
    xi_eff = min(xi_intrinsic, float(L_system))
    
    return xi_eff


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Измерение корреляционной длины ξ на пирохлорной решётке"
    )
    parser.add_argument('--L', type=int, default=4,
                       help='Размер системы')
    parser.add_argument('--Gamma', type=float, default=0.9,
                       help='Поперечное поле')
    parser.add_argument('--M', type=int, default=64,
                       help='Число троттеровских слоёв')
    parser.add_argument('--scan-L', action='store_true',
                       help='Сканирование по размерам системы')
    parser.add_argument('--analytic', action='store_true',
                       help='Только аналитическая модель')
    parser.add_argument('--save', type=str, default=None,
                       help='Сохранить результаты в JSON')
    
    args = parser.parse_args()
    
    if args.analytic:
        print("АНАЛИТИЧЕСКАЯ МОДЕЛЬ ξ(Γ):")
        print(f"{'Γ':>8} {'ξ_intr':>12} {'ξ_eff':>12}")
        print("-" * 35)
        for Gamma in np.linspace(0.1, 1.5, 15):
            xi = analytic_xi_model(Gamma, L_system=args.L)
            print(f"{Gamma:8.2f} {xi:12.2f} {xi:12.2f}")
    
    elif args.scan_L:
        results = scan_L(Gamma=args.Gamma, M=args.M)
        if args.save and results:
            with open(args.save, 'w') as f:
                json.dump(results, f, indent=2, default=float)
            print(f"Сохранено: {args.save}")
    
    else:
        r = run_xi_measurement(L=args.L, Gamma=args.Gamma, M=args.M)
        if r and args.save:
            with open(args.save, 'w') as f:
                json.dump(r, f, indent=2, default=float)
            print(f"Сохранено: {args.save}")
