"""
Квантовое Монте-Карло на пирохлорной решётке.
Suzuki-Trotter path integral + Wolff кластеры.
H = +J Σ z_i z_j (AFM, фрустрированный) − Γ Σ σ^x
"""
import numpy as np
from numba import njit
import sys, math, json
from datetime import datetime

# ============================================================
# Геометрия пирохлорной решётки (плоский массив)
# ============================================================
@njit
def build_pyrochlore(L):
    """N=4L³ спинов, до 6 соседей каждый."""
    N = 4 * L**3
    max_n = 6
    neigh = -np.ones((N, max_n), dtype=np.int32)
    deg = np.zeros(N, dtype=np.int32)
    
    def idx(x, y, z, s):
        return ((x % L) * L * L + (y % L) * L + (z % L)) * 4 + s
    
    edges = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    
    for x in range(L):
        for y in range(L):
            for z in range(L):
                for s1, s2 in edges:
                    i = idx(x,y,z,s1); j = idx(x,y,z,s2)
                    found = False
                    for k in range(deg[i]):
                        if neigh[i,k] == j: found = True; break
                    if not found and deg[i] < max_n:
                        neigh[i, deg[i]] = j; deg[i] += 1
                    found = False
                    for k in range(deg[j]):
                        if neigh[j,k] == i: found = True; break
                    if not found and deg[j] < max_n:
                        neigh[j, deg[j]] = i; deg[j] += 1
    return neigh, deg


# ============================================================
# Квантовое Монте-Карло (троттеризация)
# ============================================================
@njit
def setup_trotter(J, Gamma, M):
    """Эффективные константы после троттеризации."""
    K_spin = J / M           # пространственное AFM-взаимодействие
    x = Gamma / M
    K_tau = -0.5 * np.log(np.tanh(x)) if x > 0 else 10.0  # мнимое время (FM)
    return K_spin, K_tau


@njit
def wolff_cluster_pyro(z, neigh, deg, M, K_spin, K_tau):
    """
    Кластерное обновление Вольфа для пирохлорной решётки.
    z: [N, M] — спины (N пространственных × M троттеровских слоёв)
    """
    N = z.shape[0]
    
    # Конвертируем в ФМ переменные: s_i = (-1)^i z_i (для АФМ)
    s = np.zeros((N, M), dtype=np.float64)
    for i in range(N):
        sign = 1.0 if i % 2 == 0 else -1.0
        for tau in range(M):
            s[i, tau] = sign * z[i, tau]
    
    # Вероятности добавления связи
    p_spin = 1.0 - np.exp(-2.0 * abs(K_spin))  # пространственная (FM после s-преобразования)
    p_tau = 1.0 - np.exp(-2.0 * abs(K_tau))    # троттеровская (FM)
    
    total = N * M
    seed = np.random.randint(0, total)
    seed_tau = seed % M
    seed_i = seed // M
    
    cluster = np.zeros((N, M), dtype=np.int8)
    cluster[seed_i, seed_tau] = 1
    cluster_size = 1
    queue = [(seed_i, seed_tau)]
    head = 0
    
    while head < len(queue):
        i, tau = queue[head]; head += 1
        
        # Пространственные соседи
        for k in range(deg[i]):
            j = neigh[i, k]
            if cluster[j, tau] == 0 and s[i, tau] == s[j, tau]:
                if np.random.random() < p_spin:
                    cluster[j, tau] = 1; cluster_size += 1
                    queue.append((j, tau))
        
        # Троттеровские соседи (вдоль мнимого времени)
        tau_next = (tau + 1) % M
        if cluster[i, tau_next] == 0 and s[i, tau] == s[i, tau_next]:
            if np.random.random() < p_tau:
                cluster[i, tau_next] = 1; cluster_size += 1
                queue.append((i, tau_next))
        
        tau_prev = (tau - 1) % M
        if cluster[i, tau_prev] == 0 and s[i, tau] == s[i, tau_prev]:
            if np.random.random() < p_tau:
                cluster[i, tau_prev] = 1; cluster_size += 1
                queue.append((i, tau_prev))
    
    # Переворот кластера
    for i in range(N):
        sign = 1.0 if i % 2 == 0 else -1.0
        for tau in range(M):
            if cluster[i, tau]:
                s[i, tau] = -s[i, tau]
                z[i, tau] = sign * s[i, tau]
    
    return cluster_size


@njit
def compute_obs_qpyro(z, neigh, deg, M, K_spin, K_tau, L):
    """Наблюдаемые для квантовой пирохлорной решётки."""
    N = z.shape[0]
    
    # Энергия на узел
    E = 0.0
    for i in range(N):
        for tau in range(M):
            # Пространственные связи
            for k in range(deg[i]):
                j = neigh[i, k]
                if j > i:
                    E += K_spin * z[i, tau] * z[j, tau]
            # Троттеровские связи
            E -= K_tau * z[i, tau] * z[i, (tau+1) % M]
    E /= (N * M)
    
    # Намагниченность
    v = np.mean(z)
    
    # Staggered magnetization (АФМ параметр порядка)
    v_stag = 0.0
    for i in range(N):
        sign = 1.0 if (i % 4) < 2 else -1.0
        for tau in range(M):
            v_stag += sign * z[i, tau]
    v_stag = abs(v_stag) / (N * M)
    
    # Ice rule fraction (на каждом троттеровском слое)
    ice = 0.0
    for tau in range(M):
        ice_tau = 0
        for x in range(L):
            for y in range(L):
                for zc in range(L):
                    base = ((x*L + y)*L + zc) * 4
                    s = z[base, tau] + z[base+1, tau] + z[base+2, tau] + z[base+3, tau]
                    if abs(s) < 0.01:
                        ice_tau += 1
        ice += ice_tau
    ice /= (M * L * L * L)
    
    return E, v, v_stag, ice


# ============================================================
# Основная симуляция
# ============================================================
def run_quantum_pyro(L=2, M_trotter=16, J=1.0, Gamma=0.5, n_thermal=1000, n_samples=2000):
    """Квантовое MC на пирохлорной решётке."""
    N = 4 * L**3
    neigh, deg = build_pyrochlore(L)
    K_spin, K_tau = setup_trotter(J, Gamma, M_trotter)
    
    print(f"QMC Pyrochlore: L={L}, N={N}, M={M_trotter}, J={J}, Γ={Gamma}")
    print(f"  K_spin={K_spin:.4f}, K_tau={K_tau:.4f}, N_spins={N*M_trotter}")
    
    # Инициализация: случайная
    np.random.seed(42 + int(Gamma * 100))
    z = np.random.choice(np.array([-1.0, 1.0]), size=(N, M_trotter))
    
    # Термализация
    print(f"  Термализация ({n_thermal} кластерных шагов)...")
    for step in range(n_thermal):
        wolff_cluster_pyro(z, neigh, deg, M_trotter, K_spin, K_tau)
    
    # Измерения
    print(f"  Измерения ({n_samples} шагов)...")
    n_meas = n_samples // 5
    Es = np.zeros(n_meas); vs = np.zeros(n_meas)
    vss = np.zeros(n_meas); ices = np.zeros(n_meas)
    
    for k in range(n_meas):
        for _ in range(5):
            wolff_cluster_pyro(z, neigh, deg, M_trotter, K_spin, K_tau)
        E, v, vs_k, ice = compute_obs_qpyro(z, neigh, deg, M_trotter, K_spin, K_tau, L)
        Es[k] = E; vs[k] = v; vss[k] = vs_k; ices[k] = ice
    
    v_stag_mean = np.mean(vss)
    ice_mean = np.mean(ices)
    
    if ice_mean > 0.5 and v_stag_mean < 0.3:
        phase = "U(1)-QUANTUM-SPIN-LIQUID"
    elif v_stag_mean > 0.3:
        phase = "AFM-ordered"
    else:
        phase = "quantum-paramagnet"
    
    print(f"  E/N={np.mean(Es):.4f}±{np.std(Es):.4f}")
    print(f"  |v|={np.mean(np.abs(vs)):.4f} v_stag={v_stag_mean:.4f}±{np.std(vss):.4f}")
    print(f"  ice={ice_mean:.4f} → {phase}")
    
    return {"L":L,"N":N,"M":M_trotter,"Gamma":Gamma,
            "E":np.mean(Es),"E_std":np.std(Es),
            "v_stag":v_stag_mean,"v_stag_std":np.std(vss),
            "ice":ice_mean,"phase":phase,
            "K_spin":K_spin,"K_tau":K_tau}


# ============================================================
if __name__ == "__main__":
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    M = int(sys.argv[2]) if len(sys.argv) > 2 else 16
    
    print("="*65)
    print(f"КВАНТОВОЕ MC: ПИРОХЛОРНАЯ РЕШЁТКА")
    print(f"Поиск квантового фазового перехода Z₂→U(1)")
    print("="*65)
    
    results = []
    for Gamma in [0.1, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]:
        print(f"\n--- Γ = {Gamma:.1f} ---")
        r = run_quantum_pyro(L=L, M_trotter=M, J=1.0, Gamma=Gamma,
                            n_thermal=1000, n_samples=2000)
        results.append(r)
    
    print(f"\n{'='*65}")
    print(f"СКАНИРОВАНИЕ КВАНТОВОГО ФАЗОВОГО ПЕРЕХОДА")
    print(f"{'Γ':>6} {'v_stag':>10} {'ice':>10} {'E/N':>10} {'Фаза':>25}")
    print(f"{'─'*61}")
    for r in results:
        print(f"{r['Gamma']:6.1f} {r['v_stag']:10.4f} {r['ice']:10.4f} "
              f"{r['E']:10.4f} {r['phase']:>25}")
    
    # Сохранение
    with open(f"qmc_pyro_L{L}_M{M}.json", "w") as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\nРезультаты сохранены: qmc_pyro_L{L}_M{M}.json")
    print("КРУПНОМАСШТАБНАЯ СИМУЛЯЦИЯ ЗАПУЩЕНА.")
