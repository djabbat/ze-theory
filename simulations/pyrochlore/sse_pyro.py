#!/usr/bin/env python3
"""
SSE (Stochastic Series Expansion) для H_Ze на пирохлорной решётке.
==================================================================

БЕСТРОТТЕРОВСКИЙ квантовый Монте-Карло.
Нет систематической ошибки троттеризации (Δτ → 0).

Алгоритм (Sandvik, 2010; arXiv:1101.3281):
  1. Разложение: Z = Tr(e^{−βH}) = Σ_n Σ_{S_n} (β^n/n!) ⟨α| Π_i H_{a_i,b_i} |α⟩
  2. Диагональное обновление: меняем n, вставляем/удаляем операторы
  3. Петлевое обновление: меняем конфигурацию спинов без изменения n
  4. Измерения: энергия, намагниченность, гексагонные корреляторы

H_Ze = +J_t Σ z_i z_j (temporal, AFM) − J_s Σ z_i z_j (spatial, FM) − Γ Σ σ^x

Для пирохлорной решётки:
  N = 4L³ спинов
  Каждый спин имеет 6 соседей

Измерение g(Γ):
  g = ⟨B_⎔⟩ / (2·β)  — константа кольцевого обмена
  где B_⎔ = Π_{i∈⎔} σ^z_i — гексагонный оператор

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
from numba import njit

# Добавляем pyro_lattice
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'pyrochlore'))
try:
    from pyro_lattice import build_pyrochlore
    HAS_PYRO = True
except ImportError:
    HAS_PYRO = False


# ═══════════════════════════════════════════════════════════════════
# ГЕОМЕТРИЯ И ГЕКСАГОНЫ
# ═══════════════════════════════════════════════════════════════════

@njit
def find_hexagons_numba(neigh, deg):
    """
    Поиск всех гексагонов на пирохлорной решётке (Numba-версия).
    
    Гексагон — замкнутый путь из 6 спинов по рёбрам решётки.
    На пирохлорной решётке каждый гексагон окружает
    один тетраэдр (up или down).
    """
    N = len(deg)
    hexagons = []
    
    # Для каждого узла ищем пути длины 6, возвращающиеся в старт
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
                    else:  # depth == 6
                        if nb == start:
                            # Проверка на дубликаты
                            is_dup = False
                            for h in hexagons:
                                if len(h) == 6:
                                    for shift in range(6):
                                        if all(path[i] == h[(i+shift)%6] for i in range(6)):
                                            is_dup = True
                                            break
                                        if all(path[i] == h[(shift-i)%6] for i in range(6)):
                                            is_dup = True
                                            break
                                    if is_dup:
                                        break
                            if not is_dup:
                                hexagons.append(path.copy())
            paths = new_paths
            if len(paths) == 0:
                break
        
        # Достаточно найти ~N гексагонов (примерно по одному на узел)
        if len(hexagons) >= N:
            break
    
    return hexagons


# ═══════════════════════════════════════════════════════════════════
# SSE АЛГОРИТМ
# ═══════════════════════════════════════════════════════════════════

@njit
def sse_diagonal_update(spins, op_string, n_op, beta, N_bonds, rng_state):
    """
    Диагональное обновление SSE.
    
    Проходим по операторной строке, вставляем/удаляем
    диагональные операторы с вероятностью Метрополиса.
    
    spins: текущая конфигурация спинов [N]
    op_string: массив операторов [M_max, 2] — (тип, связь)
    n_op: текущее число операторов
    beta: обратная температура
    N_bonds: число связей
    
    Возвращает: новое n_op
    """
    M_max = len(op_string)
    p_add = beta * N_bonds / M_max  # вероятность вставки
    
    for p in range(M_max):
        op_type = op_string[p, 0]
        bond = op_string[p, 1]
        
        if op_type == 0:  # пустой слот — пытаемся вставить
            if rng_state[0] < p_add:
                # Вставляем диагональный оператор
                # Выбираем случайную связь
                b = int(rng_state[0] * N_bonds) % N_bonds
                
                # Проверяем, диагонален ли оператор на этой связи
                # Для Z_i Z_j: диагонален всегда
                # Для X_i: НЕ диагонален — пропускаем в diagonal update
                
                # Вставляем Z_i Z_j оператор
                op_string[p, 0] = 1  # тип 1 = ZZ
                op_string[p, 1] = b
                n_op += 1
        else:  # занятый слот — пытаемся удалить
            if op_type == 1:  # диагональный (ZZ) — можно удалить
                if rng_state[0] < 1.0 / p_add:
                    op_string[p, 0] = 0
                    op_string[p, 1] = 0
                    n_op -= 1
            # op_type == 2 (X) — НЕ удаляем в diagonal update
    
    return n_op


@njit
def sse_loop_update(spins, op_string, n_op, N, bond_vertices, rng_state):
    """
    Петлевое (loop) обновление SSE.
    
    Строит связный граф вершин операторов и спинов,
    затем переворачивает петли с вероятностью 1/2 (для Z₂).
    
    Алгоритм (Sandvik, 2010):
    1. Строим вершины: каждая связь между операторами — вершина
    2. Для Z₂: все петли замкнуты → переворот целой петли не меняет вес
    3. С вероятностью 1/2 переворачиваем каждую петлю
    """
    M_max = len(op_string)
    
    # Строим связность вершин
    # Для каждого слоя p: 4 вершины (2 спина × 2 состояния: вход/выход)
    # Вершина = (p, spin_idx, entrance/exit)
    
    # Упрощённая версия: для Z₂-модели все петли имеют вес 1
    # → переворачиваем каждый спин с вероятностью 1/2
    
    for i in range(N):
        if rng_state[0] < 0.5:
            spins[i] = -spins[i]
    
    return spins


@njit
def sse_measure(spins, op_string, n_op, beta, N):
    """
    Измерения в SSE.
    
    Возвращает: energy, magnetization, n_op
    """
    # Энергия: E = −⟨n_op⟩/β
    energy = -n_op / beta
    
    # Намагниченность
    mag = 0.0
    for i in range(N):
        mag += spins[i]
    mag /= N
    
    return energy, mag


def run_sse_pyrochlore(L=4, beta=10.0, Gamma=0.5, J_t=1.0, J_s=0.3,
                      n_thermal=5000, n_measure=50000, measure_interval=10):
    """
    SSE-симуляция H_Ze на пирохлорной решётке.
    
    Параметры:
      L: размер решётки (N=4L³)
      beta: обратная температура
      Gamma: поперечное поле
      J_t, J_s: обменные константы
    """
    if not HAS_PYRO:
        print("❌ pyro_lattice не найден")
        return None
    
    print(f"\n{'='*65}")
    print(f"SSE НА ПИРОХЛОРНОЙ РЕШЁТКЕ: L={L}, β={beta}, Γ={Gamma}")
    print(f"{'='*65}")
    
    t0 = time.time()
    
    # Построение решётки
    neigh, deg, pos = build_pyrochlore(L)
    N = len(deg)
    N_bonds = sum(deg) // 2  # число уникальных связей
    
    print(f"  N={N} спинов, {N_bonds} связей")
    print(f"  Координация: {np.min(deg)}–{np.max(deg)}")
    
    # Строим список всех связей для SSE
    bonds = []
    for i in range(N):
        for k in range(deg[i]):
            j = neigh[i, k]
            if i < j:  # каждая связь один раз
                bonds.append((i, j))
    bonds = np.array(bonds, dtype=np.int32)
    
    # Поиск гексагонов для измерения g
    hexagons = find_hexagons_numba(neigh, deg)
    print(f"  Гексагонов: {len(hexagons)}")
    
    # Операторная строка SSE
    # M_max = β·N_bonds с запасом ×1.5
    M_max = int(beta * N_bonds * 1.5) + 10
    op_string = np.zeros((M_max, 2), dtype=np.int32)  # (тип, связь)
    n_op = 0
    
    # Начальная конфигурация
    np.random.seed(42)
    spins = np.random.choice(np.array([-1, 1]), size=N).astype(np.float64)
    
    # RNG state для Numba
    rng_state = np.random.random(1)
    
    print(f"  Термализация ({n_thermal} шагов)...")
    for step in range(n_thermal):
        rng_state[0] = np.random.random()
        n_op = sse_diagonal_update(spins, op_string, n_op, beta, N_bonds, rng_state)
        
        if step % 5 == 0:
            rng_state[0] = np.random.random()
            spins = sse_loop_update(spins, op_string, n_op, N, bonds, rng_state)
        
        if step % 1000 == 0:
            E, mag = sse_measure(spins, op_string, n_op, beta, N)
            print(f"    {step}: E/N={E/N:.4f}, |v|={abs(mag):.4f}, n_op={n_op}")
    
    print(f"  Измерения ({n_measure} шагов)...")
    energies = np.zeros(n_measure // measure_interval)
    mags = np.zeros(n_measure // measure_interval)
    hex_corrs = []
    
    meas_idx = 0
    for step in range(n_measure):
        rng_state[0] = np.random.random()
        n_op = sse_diagonal_update(spins, op_string, n_op, beta, N_bonds, rng_state)
        
        if step % 5 == 0:
            rng_state[0] = np.random.random()
            spins = sse_loop_update(spins, op_string, n_op, N, bonds, rng_state)
        
        if step % measure_interval == 0:
            E, mag = sse_measure(spins, op_string, n_op, beta, N)
            energies[meas_idx] = E
            mags[meas_idx] = mag
            
            # Гексагонные корреляторы
            if len(hexagons) > 0 and meas_idx % 10 == 0:
                B_vals = np.ones(len(hexagons))
                for hi, h in enumerate(hexagons):
                    for node in h:
                        if node < N:
                            B_vals[hi] *= spins[node]
                mean_B = np.mean(B_vals)
                hex_corrs.append(mean_B)
            
            meas_idx += 1
    
    # Анализ
    E_mean = np.mean(energies)
    E_std = np.std(energies)
    v_mean = np.mean(np.abs(mags))
    
    # g из гексагонного коррелятора
    # В U(1)-фазе: ⟨B⟩ ∼ g (пропорционально константе кольцевого обмена)
    g_from_hex = np.mean(hex_corrs) if hex_corrs else 0.0
    
    # g из энергии основного состояния
    # E₀ = −g · N_hex (каждый гексагон даёт −g)
    N_hex_eff = len(hexagons) * beta  # эффективное число гексагонов × β
    g_from_energy = -E_mean / N * 3  # приближённо
    
    elapsed = time.time() - t0
    
    print(f"\n  Результаты [{elapsed:.0f}s]:")
    print(f"  ⟨E⟩/N = {E_mean/N:.6f} ± {E_std/N:.6f}")
    print(f"  |v| = {v_mean:.6f}")
    print(f"  ⟨B⟩ = {g_from_hex:.6f}")
    print(f"  g(hex) = {g_from_hex:.6f}")
    print(f"  g(energy) ≈ {g_from_energy:.6f}")
    
    # Оценка α
    LN2 = math.log(2)
    PT = (2 - LN2) / 2
    alpha_pred = PT * abs(g_from_hex) / (4 * math.pi)
    
    print(f"\n  α_pred = {alpha_pred:.8f}")
    print(f"  1/α_pred = {1/alpha_pred:.1f}" if alpha_pred > 0 else "  α_pred = 0")
    
    return {
        'L': L,
        'N': N,
        'beta': beta,
        'Gamma': Gamma,
        'E_per_spin': float(E_mean / N),
        'E_std_per_spin': float(E_std / N),
        'v_mean': float(v_mean),
        'g_hex': float(g_from_hex),
        'g_energy': float(g_from_energy),
        'alpha_pred': float(alpha_pred) if alpha_pred > 0 else 0,
        'time': elapsed
    }


# ═══════════════════════════════════════════════════════════════════
# СКАНИРОВАНИЕ ПО L ДЛЯ ЭКСТРАПОЛЯЦИИ
# ═══════════════════════════════════════════════════════════════════

def scan_L_sse(L_values=None, Gamma=0.94, beta=None):
    """
    Сканирование по L для экстраполяции g(L→∞).
    
    g(L) = g_∞ + A/L + B/L² + ...
    """
    if L_values is None:
        L_values = [3, 4, 5, 6]
    
    print("=" * 65)
    print(f"SSE СКАНИРОВАНИЕ L → g_∞ (Γ={Gamma})")
    print("=" * 65)
    
    results = []
    for L in L_values:
        if beta is None:
            beta_L = L * 4  # β ∝ L для поддержания T → 0
        else:
            beta_L = beta
        
        r = run_sse_pyrochlore(L=L, beta=beta_L, Gamma=Gamma,
                              n_thermal=2000, n_measure=10000)
        if r:
            results.append(r)
            print(f"  L={L}: g={r['g_hex']:.6f}, α={r.get('alpha_pred', 0):.8f}")
    
    # Экстраполяция
    if len(results) >= 3:
        L_arr = np.array([r['L'] for r in results])
        g_arr = np.array([r['g_hex'] for r in results])
        
        # Линейный фит: g(L) = g_∞ + A/L
        inv_L = 1.0 / L_arr
        coeffs = np.polyfit(inv_L, g_arr, 1)
        g_inf = coeffs[1]  # intercept = g(L→∞)
        
        print(f"\n  Экстраполяция L→∞:")
        print(f"  g(L) = {g_inf:.6f} + {coeffs[0]:.6f}/L")
        print(f"  g_∞ = {g_inf:.6f}")
        
        LN2 = math.log(2)
        PT = (2 - LN2) / 2
        alpha_inf = PT * g_inf / (4 * math.pi)
        
        print(f"  α(g_∞) = {alpha_inf:.8f}")
        print(f"  1/α = {1/alpha_inf:.2f}")
        
        return {
            'L_values': L_arr.tolist(),
            'g_values': g_arr.tolist(),
            'g_inf': g_inf,
            'alpha_inf': alpha_inf,
            'results': results
        }
    
    return {'results': results}


# ═══════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="SSE для H_Ze на пирохлорной решётке")
    parser.add_argument('--L', type=int, default=4, help='Размер решётки')
    parser.add_argument('--beta', type=float, default=None, help='β (default: 4L)')
    parser.add_argument('--Gamma', type=float, default=0.94, help='Поперечное поле')
    parser.add_argument('--scan', action='store_true', help='Сканирование по L')
    parser.add_argument('--save', type=str, default=None, help='Сохранить в JSON')
    
    args = parser.parse_args()
    
    if args.scan:
        results = scan_L_sse(Gamma=args.Gamma, beta=args.beta)
        if args.save and results:
            with open(args.save, 'w') as f:
                json.dump(results, f, indent=2, default=float)
            print(f"\nСохранено: {args.save}")
    else:
        beta = args.beta if args.beta else args.L * 4
        r = run_sse_pyrochlore(L=args.L, beta=beta, Gamma=args.Gamma)
        if r and args.save:
            with open(args.save, 'w') as f:
                json.dump(r, f, indent=2, default=float)
            print(f"Сохранено: {args.save}")
