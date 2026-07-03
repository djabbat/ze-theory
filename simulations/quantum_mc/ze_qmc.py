"""
Квантовое Монте-Карло для Ze-гамильтониана с поперечным полем.

H = -J_t Σ z_t z_{t+1} - J_s Σ_{⟨x,y⟩} z_x z_y - Γ Σ σ^x

Метод: Suzuki-Trotter path integral → отображение в классическую
модель Изинга в (d+1) измерениях (добавляется мнимое время).

J_τ = -½ ln tanh(Γ / M_trotter)

Алгоритм: Метрополис + кластерные обновления (Wolff) для
преодоления критического замедления вблизи фазового перехода.
"""

import numpy as np
from numba import njit
from pathlib import Path
import json
from datetime import datetime


DEFAULT_PARAMS = {
    "L": 16,                # пространственно-временной размер (1+1d цепочка)
    "M_trotter": 32,        # число троттеровских слоёв (мнимое время)
    "J_t": 1.0,             # антиферромагнитная связь
    "Gamma": 1.0,           # поперечное поле (квантовые флуктуации)
    "h": 0.0,               # продольное поле
    "n_thermal": 5000,      # термализация
    "n_samples": 20000,     # измерения
    "sample_interval": 10,  # интервал
    "seed": 42,
}


@njit
def setup_trotter(J_t, Gamma, M, h):
    """Вычисление эффективных констант связей."""
    # Взаимодействие вдоль реального времени (антиферромагнитное)
    K_t = J_t / M  # делим на M из-за троттеризации
    
    # Взаимодействие вдоль мнимого времени (ферромагнитное)
    # J_τ = -½ ln tanh(Γ/M)
    x = Gamma / M
    if x > 0:
        K_tau = -0.5 * np.log(np.tanh(x))
    else:
        K_tau = 10.0  # большой, но конечный (Γ=0 → нет квантовых флуктуаций)
    
    # Поле (делим на M)
    K_h = h / M
    
    return K_t, K_tau, K_h


@njit
def wolff_cluster_update(z, L, M, K_t, K_tau, K_h):
    """Кластерное обновление Вольфа для антиферромагнетика.
    
    Адаптировано для 2D решётки (L реального времени × M мнимого).
    Антиферромагнитное взаимодействие вдоль реального времени (K_t > 0)
    означает, что спины на соседних t-узлах стремятся быть противоположными.
    
    Для кластерного алгоритма переходим к переменным s_i = (-1)^i z_i
    чтобы превратить АФМ взаимодействие в ФМ.
    """
    # Конвертируем в ферромагнитные переменные
    # s_{t,τ} = (-1)^t * z_{t,τ}
    s = np.zeros((L, M), dtype=np.float64)
    for t in range(L):
        sign = 1.0 if t % 2 == 0 else -1.0
        for tau in range(M):
            s[t, tau] = sign * z[t, tau]
    
    total = L * M
    # Выбираем случайный затравочный узел
    seed_idx = np.random.randint(0, total)
    seed_tau = seed_idx % M
    tmp = seed_idx // M
    seed_t = tmp % L
    
    # Вероятность добавления связи (ферромагнитная, K > 0)
    p_add_t = 1.0 - np.exp(-2.0 * abs(K_t))  # вдоль реального времени
    p_add_tau = 1.0 - np.exp(-2.0 * abs(K_tau))  # вдоль мнимого времени
    
    # Построение кластера (BFS)
    cluster = np.zeros((L, M), dtype=np.int8)
    cluster[seed_t, seed_tau] = 1
    cluster_size = 1
    
    queue = [(seed_t, seed_tau)]
    queue_head = 0
    
    while queue_head < len(queue):
        t, tau = queue[queue_head]
        queue_head += 1
        
        # Сосед вдоль реального времени (вперёд)
        t_next = (t + 1) % L
        if cluster[t_next, tau] == 0 and s[t, tau] == s[t_next, tau]:
            if np.random.random() < p_add_t:
                cluster[t_next, tau] = 1
                cluster_size += 1
                queue.append((t_next, tau))
        
        # Сосед вдоль реального времени (назад)
        t_prev = (t - 1) % L
        if cluster[t_prev, tau] == 0 and s[t, tau] == s[t_prev, tau]:
            if np.random.random() < p_add_t:
                cluster[t_prev, tau] = 1
                cluster_size += 1
                queue.append((t_prev, tau))
        
        # Сосед вдоль мнимого времени (вперёд)
        tau_next = (tau + 1) % M
        if cluster[t, tau_next] == 0 and s[t, tau] == s[t, tau_next]:
            if np.random.random() < p_add_tau:
                cluster[t, tau_next] = 1
                cluster_size += 1
                queue.append((t, tau_next))
        
        # Сосед вдоль мнимого времени (назад)
        tau_prev = (tau - 1) % M
        if cluster[t, tau_prev] == 0 and s[t, tau] == s[t, tau_prev]:
            if np.random.random() < p_add_tau:
                cluster[t, tau_prev] = 1
                cluster_size += 1
                queue.append((t, tau_prev))
    
    # Переворачиваем кластер
    for t in range(L):
        sign = 1.0 if t % 2 == 0 else -1.0
        for tau in range(M):
            if cluster[t, tau]:
                s[t, tau] = -s[t, tau]
                z[t, tau] = sign * s[t, tau]
    
    return cluster_size


@njit
def metropolis_sweep(z, L, M, K_t, K_tau, K_h, n_sweeps=1):
    """Одночастичные обновления Метрополиса (дополнительно к Вольфу)."""
    total = L * M
    accepted = 0
    attempts = 0
    
    for _ in range(n_sweeps):
        for _ in range(total):
            idx = np.random.randint(0, total)
            tau = idx % M
            tmp = idx // M
            t = tmp % L
            
            # Локальная энергия
            e_old = -K_h * z[t, tau]
            # реальное время (антиферромагнитное)
            e_old += -K_t * z[t, tau] * z[(t-1)%L, tau]
            e_old += -K_t * z[t, tau] * z[(t+1)%L, tau]
            # мнимое время (ферромагнитное)
            e_old += -K_tau * z[t, tau] * z[t, (tau-1)%M]
            e_old += -K_tau * z[t, tau] * z[t, (tau+1)%M]
            
            z[t, tau] = -z[t, tau]
            
            e_new = -K_h * z[t, tau]
            e_new += -K_t * z[t, tau] * z[(t-1)%L, tau]
            e_new += -K_t * z[t, tau] * z[(t+1)%L, tau]
            e_new += -K_tau * z[t, tau] * z[t, (tau-1)%M]
            e_new += -K_tau * z[t, tau] * z[t, (tau+1)%M]
            
            delta = e_new - e_old
            attempts += 1
            
            if delta <= 0 or np.random.random() < np.exp(-delta):
                accepted += 1
            else:
                z[t, tau] = -z[t, tau]
    
    return accepted / max(attempts, 1)


@njit
def compute_observables(z, L, M, K_t, K_tau):
    """Вычисление физических наблюдаемых.
    
    Возвращает: энергию, намагниченность |v|,
    v_staggered (параметр порядка АФМ),
    корреляционные функции.
    """
    # Энергия на узел (в единицах исходного H)
    E = 0.0
    for t in range(L):
        for tau in range(M):
            E += -K_t * z[t, tau] * z[(t+1)%L, tau]
            E += -K_tau * z[t, tau] * z[t, (tau+1)%M]
    E /= (L * M)
    
    # Намагниченность vz = ⟨z⟩
    vz = np.mean(z)
    
    # Стаггерированная намагниченность (АФМ параметр порядка)
    # v_stag = ⟨(-1)^t z_t⟩
    v_stag = 0.0
    for t in range(L):
        sign = 1.0 if t % 2 == 0 else -1.0
        for tau in range(M):
            v_stag += sign * z[t, tau]
    v_stag /= (L * M)
    
    # Корреляционная функция вдоль реального времени (АФМ)
    max_dist = L // 2
    corr_t = np.zeros(max_dist + 1)
    for t in range(L):
        for tau in range(M):
            for d in range(max_dist + 1):
                t2 = (t + d) % L
                phase = (-1) ** d  # антиферромагнитная фаза
                corr_t[d] += phase * z[t, tau] * z[t2, tau]
    corr_t /= (L * M * (max_dist + 1) / (max_dist + 1))  # нормировка
    # правильная нормировка
    norm = L * M
    corr_t /= norm
    # пересчёт — каждое расстояние считается L*M раз
    corr_t_corrected = np.zeros(max_dist + 1)
    for d in range(max_dist + 1):
        total = 0.0
        count = 0
        for t in range(L):
            for tau in range(M):
                t2 = (t + d) % L
                phase = (-1) ** d
                total += phase * z[t, tau] * z[t2, tau]
                count += 1
        corr_t_corrected[d] = total / count
    corr_t = corr_t_corrected
    
    # Магнитная восприимчивость (флуктуации намагниченности)
    # χ = β (⟨v²⟩ - ⟨v⟩²) * N
    # но в троттеровской формулировке β эффективно = M
    # χ_scaled = M * (⟨v_stag²⟩ - ⟨v_stag⟩²) * L
    
    return E, vz, v_stag, corr_t


class QuantumZeSimulation:
    """Квантовое Монте-Карло для Ze-модели."""
    
    def __init__(self, **kwargs):
        self.params = DEFAULT_PARAMS.copy()
        self.params.update(kwargs)
        p = self.params
        self.K_t, self.K_tau, self.K_h = setup_trotter(
            p["J_t"], p["Gamma"], p["M_trotter"], p["h"]
        )
        
    def run(self):
        p = self.params
        L, M = p["L"], p["M_trotter"]
        
        print(f"Квантовое MC: L={L}, M_trotter={M}, J_t={p['J_t']}, Γ={p['Gamma']}")
        print(f"  Эффективные связи: K_t={self.K_t:.4f}, K_tau={self.K_tau:.4f}")
        
        np.random.seed(p["seed"])
        
        # Инициализация: staggered order (антиферромагнетик)
        z = np.ones((L, M), dtype=np.float64)
        for t in range(L):
            if t % 2 == 1:
                z[t, :] = -1.0
        
        # Термализация
        print(f"  Термализация: {p['n_thermal']} шагов...")
        for step in range(p['n_thermal']):
            if step % 10 == 0:
                wolff_cluster_update(z, L, M, self.K_t, self.K_tau, self.K_h)
            else:
                metropolis_sweep(z, L, M, self.K_t, self.K_tau, self.K_h, n_sweeps=1)
        
        # Измерения
        print(f"  Измерения: {p['n_samples']} шагов...")
        n_meas = p['n_samples'] // p['sample_interval']
        energies = np.zeros(n_meas)
        magnetizations = np.zeros(n_meas)
        stag_magnetizations = np.zeros(n_meas)
        
        meas_idx = 0
        for step in range(p['n_samples']):
            if step % 5 == 0:
                wolff_cluster_update(z, L, M, self.K_t, self.K_tau, self.K_h)
            metropolis_sweep(z, L, M, self.K_t, self.K_tau, self.K_h, n_sweeps=1)
            
            if step % p['sample_interval'] == 0:
                E, vz, v_stag, _ = compute_observables(z, L, M, self.K_t, self.K_tau)
                energies[meas_idx] = E
                magnetizations[meas_idx] = abs(vz)
                stag_magnetizations[meas_idx] = abs(v_stag)
                meas_idx += 1
        
        # Корреляции (на финальной конфигурации)
        _, _, _, corr_t = compute_observables(z, L, M, self.K_t, self.K_tau)
        
        # Корреляционная длина
        xi = self._correlation_length(corr_t)
        
        v_star = 1.0 - np.log(2.0)
        
        results = {
            "params": p,
            "K_t": self.K_t,
            "K_tau": self.K_tau,
            "mean_energy": float(np.mean(energies)),
            "std_energy": float(np.std(energies)),
            "mean_v": float(np.mean(magnetizations)),
            "std_v": float(np.std(magnetizations)),
            "mean_v_stag": float(np.mean(stag_magnetizations)),
            "std_v_stag": float(np.std(stag_magnetizations)),
            "xi_t": xi,
            "v_star_theory": v_star,
            "delta_v_star": float(abs(np.mean(magnetizations) - v_star)),
            "corr_t": corr_t.tolist(),
        }
        
        print(f"\n  Результаты:")
        print(f"  ⟨E⟩ = {results['mean_energy']:.4f} ± {results['std_energy']:.4f}")
        print(f"  |v| = {results['mean_v']:.4f} ± {results['std_v']:.4f}")
        print(f"  |v_stag| = {results['mean_v_stag']:.4f} ± {results['std_v_stag']:.4f}")
        print(f"  ξ_t = {xi:.2f}")
        print(f"  v* = {v_star:.4f}, Δv = {results['delta_v_star']:.4f}")
        
        # Интерпретация фазы
        if results['mean_v_stag'] > 0.3:
            print(f"  → ФАЗА: антиферромагнитный порядок (конфайнмент)")
        elif results['mean_v'] < 0.15:
            print(f"  → ФАЗА: парамагнетик / квантовый беспорядок")
        else:
            print(f"  → ФАЗА: промежуточная / критическая область")
        
        self.results = results
        self.final_config = z
        return results
    
    @staticmethod
    def _correlation_length(corr):
        c = np.abs(corr[1:])
        mask = c > 0.02
        if not mask.any():
            return 0.5
        d = np.arange(1, len(corr))[mask]
        v = c[mask]
        if len(v) < 2:
            return 0.5
        coeffs = np.polyfit(d, np.log(v + 1e-15), 1)
        return float(-1.0 / coeffs[0]) if coeffs[0] < 0 else float('inf')
    
    def save_results(self, filepath):
        output = {"timestamp": datetime.now().isoformat(), **self.results}
        with open(filepath, 'w') as f:
            json.dump(output, f, indent=2, default=lambda x: float(x) if isinstance(x, (np.floating,)) else x)
        print(f"  Сохранено: {filepath}")


# ============================================================
# Сканирование по Γ (поиск квантового фазового перехода)
# ============================================================
if __name__ == "__main__":
    import sys
    
    if "--quick" in sys.argv:
        sim = QuantumZeSimulation(L=8, M_trotter=16, Gamma=1.0, J_t=1.0,
                                   n_thermal=2000, n_samples=5000)
        sim.run()
        sim.save_results("qmc_results_quick.json")
    
    elif "--scan" in sys.argv:
        print("Сканирование квантового фазового перехода по Γ...\n")
        print(f"{'Γ':>8} {'|v|':>10} {'|v_stag|':>10} {'E':>10} {'ξ_t':>8} {'Фаза':>20}")
        print("-"*70)
        
        for Gamma in [0.2, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0]:
            sim = QuantumZeSimulation(L=12, M_trotter=24, Gamma=Gamma, J_t=1.0,
                                       n_thermal=3000, n_samples=10000)
            r = sim.run()
            
            if r['mean_v_stag'] > 0.3:
                phase = "АФМ (конфайнмент)"
            elif r['mean_v'] < 0.15:
                phase = "квантовый парамагнетик"
            else:
                phase = "критическая"
            
            print(f"{Gamma:8.2f} {r['mean_v']:10.4f} {r['mean_v_stag']:10.4f} "
                  f"{r['mean_energy']:10.4f} {r['xi_t']:8.2f} {phase:>20}")
    
    else:
        sim = QuantumZeSimulation(L=16, M_trotter=32, Gamma=1.0, J_t=1.0)
        sim.run()
        sim.save_results("qmc_results.json")
