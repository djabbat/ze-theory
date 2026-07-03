"""
Классическое Монте-Карло для Ze-гамильтониана H_Ze.

H_Ze = -J_t Σ_{x,t} z_{x,t} z_{x,t+1} - J_s Σ_{⟨x,y⟩,t} z_{x,t} z_{y,t} - h Σ_{x,t} z_{x,t}

z ∈ {+1 (T), -1 (S)}
J_t > 0 — антиферромагнетик вдоль времени (антипараллелизм S=−T)
J_s > 0 — ферромагнетик вдоль пространства (локальность)
h — внешнее поле (асимметрия T/S)

Алгоритм: Метрополис
Граничные условия: периодические
"""

import numpy as np
from numba import njit
from pathlib import Path
import json
from datetime import datetime


# ============================================================
# Параметры по умолчанию
# ============================================================
DEFAULT_PARAMS = {
    "Lx": 8,               # пространственный размер x
    "Ly": 8,               # пространственный размер y
    "Lt": 16,              # временной размер
    "J_t": 1.0,            # антиферромагнитная связь вдоль t (>0)
    "J_s": 0.3,            # ферромагнитная связь вдоль x,y (>0)
    "h": 0.0,              # внешнее поле
    "T": 2.0,              # температура (β = 1/T)
    "n_thermal": 5000,     # шагов термализации
    "n_samples": 20000,    # шагов измерения
    "sample_interval": 5,  # интервал между измерениями
    "seed": 42,            # зерно ГПСЧ
}


# ============================================================
# Ядро Монте-Карло (Numba-оптимизация)
# ============================================================
@njit
def energy_local(z, x, y, t, Lx, Ly, Lt, J_t, J_s, h):
    """Энергия одного спина и его связей."""
    Lx_val, Ly_val, Lt_val = Lx, Ly, Lt
    e = -h * z[x, y, t]
    # временные связи (антиферромагнитные — поощряют чередование)
    e += -J_t * z[x, y, t] * z[x, y, (t - 1) % Lt_val]
    e += -J_t * z[x, y, t] * z[x, y, (t + 1) % Lt_val]
    # пространственные связи x (ферромагнитные — поощряют одинаковость)
    e += -J_s * z[x, y, t] * z[(x - 1) % Lx_val, y, t]
    e += -J_s * z[x, y, t] * z[(x + 1) % Lx_val, y, t]
    # пространственные связи y
    e += -J_s * z[x, y, t] * z[x, (y - 1) % Ly_val, t]
    e += -J_s * z[x, y, t] * z[x, (y + 1) % Ly_val, t]
    return e


@njit
def total_energy(z, Lx, Ly, Lt, J_t, J_s, h):
    """Полная энергия конфигурации (с учётом двойного счёта связей / 2)."""
    E = 0.0
    for x in range(Lx):
        for y in range(Ly):
            for t in range(Lt):
                E += -h * z[x, y, t]
                E += -J_t * z[x, y, t] * z[x, y, (t + 1) % Lt]
                E += -J_s * z[x, y, t] * z[(x + 1) % Lx, y, t]
                E += -J_s * z[x, y, t] * z[x, (y + 1) % Ly, t]
    return E


@njit
def magnetization(z):
    """Ze-скорость v = ⟨z⟩ = (N_T - N_S)/(N_T + N_S)."""
    return np.mean(z)


@njit
def run_metropolis(z, Lx, Ly, Lt, J_t, J_s, h, beta, n_steps, sample_interval):
    """Метрополис-семплирование."""
    n_measurements = n_steps // sample_interval
    energies = np.zeros(n_measurements)
    magnetizations = np.zeros(n_measurements)
    
    total = Lx * Ly * Lt
    accepted = 0
    meas_idx = 0
    
    for step in range(n_steps):
        # случайный узел
        idx = np.random.randint(0, total)
        t = idx % Lt
        idx //= Lt
        y = idx % Ly
        x = idx // Ly
        
        # энергия до переворота
        e_old = energy_local(z, x, y, t, Lx, Ly, Lt, J_t, J_s, h)
        
        # переворот спина
        z[x, y, t] = -z[x, y, t]
        
        # энергия после переворота
        e_new = energy_local(z, x, y, t, Lx, Ly, Lt, J_t, J_s, h)
        
        delta_e = e_new - e_old
        
        # Метрополис
        if delta_e <= 0 or np.random.random() < np.exp(-beta * delta_e):
            accepted += 1
        else:
            z[x, y, t] = -z[x, y, t]  # отвергаем — возвращаем
        
        # измерение
        if step % sample_interval == 0:
            energies[meas_idx] = total_energy(z, Lx, Ly, Lt, J_t, J_s, h)
            magnetizations[meas_idx] = magnetization(z)
            meas_idx += 1
    
    return energies, magnetizations, accepted / n_steps


@njit
def compute_correlations(z, max_dist, Lx, Ly, Lt, J_t_sign=-1):
    """Корреляции вдоль времени и пространства.
    J_t_sign = -1 означает, что мы измеряем корреляции
    с чередующимся знаком вдоль времени (антиферромагнитный порядок).
    """
    corr_t = np.zeros(max_dist + 1)
    corr_x = np.zeros(max_dist + 1)
    count_t = np.zeros(max_dist + 1)
    count_x = np.zeros(max_dist + 1)
    
    for x in range(Lx):
        for y in range(Ly):
            for t in range(Lt):
                for d in range(max_dist + 1):
                    # временная корреляция
                    t2 = (t + d) % Lt
                    phase = J_t_sign ** d  # чередование знака для АФМ
                    corr_t[d] += phase * z[x, y, t] * z[x, y, t2]
                    count_t[d] += 1
                    
                    # пространственная корреляция (вдоль x)
                    x2 = (x + d) % Lx
                    corr_x[d] += z[x, y, t] * z[x2, y, t]
                    count_x[d] += 1
    
    corr_t /= count_t
    corr_x /= count_x
    return corr_t, corr_x


# ============================================================
# Основная симуляция
# ============================================================
class ZeSimulation:
    """Классическое Монте-Карло для Ze-модели."""
    
    def __init__(self, **kwargs):
        self.params = DEFAULT_PARAMS.copy()
        self.params.update(kwargs)
        self._validate()
        
    def _validate(self):
        p = self.params
        assert p["Lx"] > 0 and p["Ly"] > 0 and p["Lt"] > 0
        assert p["J_t"] > 0, "J_t должно быть > 0 (антиферромагнетик)"
        assert p["J_s"] > 0, "J_s должно быть > 0 (ферромагнетик)"
        assert p["T"] > 0
        
    @property
    def beta(self):
        return 1.0 / self.params["T"]
    
    def run(self):
        p = self.params
        Lx, Ly, Lt = p["Lx"], p["Ly"], p["Lt"]
        
        # инициализация (случайная)
        np.random.seed(p["seed"])
        z = np.random.choice([-1, 1], size=(Lx, Ly, Lt)).astype(np.float64)
        
        # термализация
        print(f"Термализация: {p['n_thermal']} шагов...")
        _, _, acc = run_metropolis(
            z, Lx, Ly, Lt, p["J_t"], p["J_s"], p["h"],
            self.beta, p["n_thermal"], p["n_thermal"]
        )
        print(f"  acceptance rate: {acc:.3f}")
        
        # измерения
        print(f"Измерения: {p['n_samples']} шагов...")
        energies, magnetizations, acc = run_metropolis(
            z, Lx, Ly, Lt, p["J_t"], p["J_s"], p["h"],
            self.beta, p["n_samples"], p["sample_interval"]
        )
        print(f"  acceptance rate: {acc:.3f}")
        
        # корреляции
        max_dist = min(Lt, Lx, Ly) // 2
        corr_t, corr_x = compute_correlations(
            z, max_dist, Lx, Ly, Lt,
            J_t_sign=-1  # антиферромагнитное упорядочение вдоль t
        )
        
        # результаты
        results = {
            "params": p,
            "mean_energy": float(np.mean(energies)),
            "std_energy": float(np.std(energies)),
            "mean_magnetization": float(np.mean(magnetizations)),
            "std_magnetization": float(np.std(magnetizations)),
            "magnetization_abs": float(np.mean(np.abs(magnetizations))),
            "corr_t": corr_t.tolist(),
            "corr_x": corr_x.tolist(),
            "energies": energies.tolist(),
            "magnetizations": magnetizations.tolist(),
        }
        
        # Ze-параметры
        v_mean = results["mean_magnetization"]
        v_abs = results["magnetization_abs"]
        results["ze_velocity_v"] = v_mean
        results["ze_velocity_abs"] = v_abs
        
        # оценка корреляционной длины (экспоненциальный фит)
        results["xi_t"] = self._correlation_length(corr_t)
        results["xi_x"] = self._correlation_length(corr_x)
        
        print(f"\nРезультаты:")
        print(f"  ⟨E⟩/N = {results['mean_energy'] / (Lx*Ly*Lt):.4f}")
        print(f"  v = ⟨z⟩ = {v_mean:.4f} ± {results['std_magnetization']:.4f}")
        print(f"  |v| = {v_abs:.4f}")
        print(f"  ξ_t = {results['xi_t']:.2f}, ξ_x = {results['xi_x']:.2f}")
        
        # сравнение с v* = 1 - ln 2
        v_star = 1.0 - np.log(2.0)
        print(f"  v* = {v_star:.4f}  (теоретическое)")
        print(f"  |v - v*| = {abs(v_abs - v_star):.4f}")
        
        self.results = results
        self.final_config = z
        return results
    
    @staticmethod
    def _correlation_length(corr):
        """Оценка корреляционной длины по экспоненциальному спаду."""
        # берём абсолютные значения и убираем нулевое расстояние
        c = np.abs(corr[1:])
        # находим, где корреляция падает ниже 0.1
        mask = c > 0.05
        if not mask.any():
            return 1.0
        distances = np.arange(1, len(corr))[mask]
        values = c[mask]
        if len(values) < 2:
            return 1.0
        # линейный фит log(c) ~ -d/xi
        coeffs = np.polyfit(distances, np.log(values), 1)
        xi = -1.0 / coeffs[0] if coeffs[0] < 0 else float('inf')
        return xi
    
    def save_results(self, filepath):
        """Сохранить результаты в JSON."""
        output = {
            "timestamp": datetime.now().isoformat(),
            **self.results,
        }
        # конвертируем numpy-типы
        def convert(obj):
            if isinstance(obj, (np.integer,)):
                return int(obj)
            if isinstance(obj, (np.floating,)):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return obj
        
        with open(filepath, 'w') as f:
            json.dump(output, f, indent=2, default=convert)
        print(f"\nРезультаты сохранены в {filepath}")


# ============================================================
# Сканирование по температуре (поиск фазового перехода)
# ============================================================
def scan_temperature(T_range, **kwargs):
    """Сканирование по температуре для поиска критической точки."""
    results = []
    for T in T_range:
        print(f"\n{'='*50}")
        print(f"T = {T:.3f}  (β = {1/T:.3f})")
        print(f"{'='*50}")
        sim = ZeSimulation(T=T, **kwargs)
        r = sim.run()
        results.append(r)
    return results


# ============================================================
# Запуск
# ============================================================
if __name__ == "__main__":
    import sys
    
    # Быстрый тест с малыми размерами
    if "--quick" in sys.argv:
        params = {"Lx": 4, "Ly": 4, "Lt": 8, "n_thermal": 1000, "n_samples": 5000}
        sim = ZeSimulation(**params)
        results = sim.run()
        sim.save_results("results_quick.json")
    
    # Сканирование по J_s (пространственная связь)
    elif "--scan-js" in sys.argv:
        Js_values = np.linspace(0.1, 1.0, 10)
        params = {"Lx": 6, "Ly": 6, "Lt": 12, "n_thermal": 3000, "n_samples": 10000}
        all_results = []
        for Js in Js_values:
            sim = ZeSimulation(J_s=Js, **params)
            r = sim.run()
            all_results.append(r)
        
        # сводка
        print("\n" + "="*60)
        print("Сводка сканирования по J_s:")
        print(f"{'J_s':>8} {'⟨E⟩/N':>10} {'|v|':>8} {'ξ_t':>8} {'ξ_x':>8}")
        for r in all_results:
            N = r["params"]["Lx"] * r["params"]["Ly"] * r["params"]["Lt"]
            print(f"{r['params']['J_s']:8.3f} {r['mean_energy']/N:10.4f} "
                  f"{r['magnetization_abs']:8.4f} {r['xi_t']:8.2f} {r['xi_x']:8.2f}")
    
    # Полноразмерное сканирование по температуре
    else:
        print("Запуск полной симуляции...")
        sim = ZeSimulation(Lx=8, Ly=8, Lt=16, J_s=0.3, T=2.0,
                          n_thermal=5000, n_samples=20000)
        results = sim.run()
        sim.save_results("results_full.json")
