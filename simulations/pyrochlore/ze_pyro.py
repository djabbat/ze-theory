"""
Ze на пирохлорной решётке: поиск U(1)-фазы.
Numba-совместимая версия.
"""
import numpy as np
from numba import njit
import sys, math
from datetime import datetime

# ============================================================
@njit
def build_pyrochlore_flat(L):
    """
    Пирохлорная решётка: N=4L³ спинов, 6 соседей на спин.
    Возвращает: neighbors[N,6] — индексы соседей.
    """
    N = 4 * L**3
    max_neigh = 6
    neigh = -np.ones((N, max_neigh), dtype=np.int32)
    deg = np.zeros(N, dtype=np.int32)
    
    def idx(x, y, z, s):
        return ((x % L) * L * L + (y % L) * L + (z % L)) * 4 + s
    
    tetra_edges = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    
    for x in range(L):
        for y in range(L):
            for z in range(L):
                for s1, s2 in tetra_edges:
                    i = idx(x,y,z,s1)
                    j = idx(x,y,z,s2)
                    # Добавляем если ещё нет
                    found = False
                    for k in range(deg[i]):
                        if neigh[i,k] == j:
                            found = True
                            break
                    if not found and deg[i] < max_neigh:
                        neigh[i, deg[i]] = j
                        deg[i] += 1
                    found = False
                    for k in range(deg[j]):
                        if neigh[j,k] == i:
                            found = True
                            break
                    if not found and deg[j] < max_neigh:
                        neigh[j, deg[j]] = i
                        deg[j] += 1
    
    return neigh, deg


@njit
def energy_pyro(z, neigh, deg, J):
    E = 0.0
    N = len(z)
    for i in range(N):
        for k in range(deg[i]):
            j = neigh[i,k]
            if j > i:
                E += J * z[i] * z[j]
    return E


@njit
def mc_step_pyro(z, neigh, deg, J, T, n_steps):
    N = len(z)
    beta = 1.0 / T
    accepted = 0
    
    for step in range(n_steps):
        i = np.random.randint(0, N)
        e_old = 0.0
        for k in range(deg[i]):
            e_old += J * z[i] * z[neigh[i,k]]
        z[i] = -z[i]
        e_new = 0.0
        for k in range(deg[i]):
            e_new += J * z[i] * z[neigh[i,k]]
        delta = e_new - e_old
        if delta <= 0 or np.random.random() < np.exp(-beta * delta):
            accepted += 1
        else:
            z[i] = -z[i]
    return accepted / (N * n_steps)


@njit
def compute_obs_pyro(z, neigh, deg, J, L):
    N = len(z)
    m = 0.0
    for i in range(N):
        m += z[i]
    m /= N
    
    E = 0.0
    for i in range(N):
        for k in range(deg[i]):
            j = neigh[i,k]
            if j > i:
                E += J * z[i] * z[j]
    E /= N
    
    # Staggered magnetization
    v_stag = 0.0
    for i in range(N):
        s = i % 4
        sign = 1.0 if s < 2 else -1.0
        v_stag += sign * z[i]
    v_stag = abs(v_stag) / N
    
    # Ice rule fraction
    ice_count = 0
    for x in range(L):
        for y in range(L):
            for zc in range(L):
                base = ((x*L + y)*L + zc) * 4
                s = z[base] + z[base+1] + z[base+2] + z[base+3]
                if abs(s) < 0.01:
                    ice_count += 1
    ice_frac = ice_count / (L*L*L)
    
    return E, m, v_stag, ice_frac


def run_pyro(L=3, J=1.0, T=0.5, n_thermal=3000, n_samples=5000):
    N = 4 * L**3
    print(f"Pyrochlore L={L}, N={N}, J={J}, T={T}")
    
    neigh, deg = build_pyrochlore_flat(L)
    print(f"  Координация: mean={np.mean(deg):.1f}, min={np.min(deg)}, max={np.max(deg)}")
    
    np.random.seed(42)
    z = np.random.choice(np.array([-1.0, 1.0]), size=N)
    
    # Термализация
    acc = mc_step_pyro(z, neigh, deg, J, T, n_thermal)
    print(f"  Термализация: acc={acc:.3f}")
    
    # Измерения
    n_meas = n_samples // 10
    Es = np.zeros(n_meas); ms = np.zeros(n_meas)
    vs = np.zeros(n_meas); ices = np.zeros(n_meas)
    
    for k in range(n_meas):
        mc_step_pyro(z, neigh, deg, J, T, 10)
        E, m, vs_k, ice = compute_obs_pyro(z, neigh, deg, J, L)
        Es[k] = E; ms[k] = m; vs[k] = vs_k; ices[k] = ice
    
    v_stag_mean = np.mean(vs)
    ice_mean = np.mean(ices)
    
    if ice_mean > 0.5 and v_stag_mean < 0.3:
        phase = "U(1)-SPIN-LIQUID"
    elif v_stag_mean > 0.3:
        phase = "AFM"
    else:
        phase = "paramagnet"
    
    print(f"  E/N={np.mean(Es):.4f}±{np.std(Es):.4f} |m|={np.mean(np.abs(ms)):.4f} v_stag={v_stag_mean:.4f} ice={ice_mean:.4f} → {phase}")
    
    return {"L":L,"N":N,"J":J,"T":T,"E":np.mean(Es),"E_std":np.std(Es),
            "v_stag":v_stag_mean,"v_stag_std":np.std(vs),
            "ice":ice_mean,"phase":phase}


if __name__ == "__main__":
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    print("="*60)
    print("Ze на пирохлорной решётке — поиск U(1)-фазы")
    print("="*60)
    
    for T in [0.2, 0.5, 1.0, 2.0, 5.0]:
        print(f"\n--- T = {T:.1f} ---")
        run_pyro(L=L, J=1.0, T=T, n_thermal=2000, n_samples=5000)
    
    print("\nU(1)-фаза: ice→1, v_stag→0 (деконфайнмент)")
    print("АФМ: v_stag→1 (конфайнмент)")
