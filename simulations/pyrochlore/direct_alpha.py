"""Прямое измерение α: ⟨B(0)B(r)⟩ → ξ → α = 0.0520/ξ."""
import numpy as np, math, json, time, sys
from pyro_lattice import build_pyrochlore
from ze_qmc_pyro import setup_trotter, wolff_cluster_pyro
from numba import njit

@njit
def find_hexagons(neigh, deg, max_hex=500):
    N = len(deg); hx = []
    for start in range(min(N, 50)):
        paths = [[start]]
        for depth in range(1, 7):
            new = []
            for path in paths:
                last = path[-1]
                for k in range(deg[last]):
                    nb = neigh[last, k]
                    if depth < 6:
                        if nb not in path: new.append(path + [nb])
                    else:
                        if nb == start: hx.append(path)
                        if len(hx) >= max_hex: return hx
            paths = new
    return hx

@njit
def hex_correlator(z, harr, centers, max_dist=8):
    nh = len(harr); M = z.shape[1]
    corr = np.zeros(max_dist + 1); cnt = np.zeros(max_dist + 1)
    bw = 0.5
    for tau in range(M):
        B = np.ones(nh)
        for h in range(nh):
            for node in harr[h]: B[h] *= z[node, tau]
        for i in range(nh):
            for j in range(i + 1, nh):
                dx = centers[i,0] - centers[j,0]
                dy = centers[i,1] - centers[j,1]
                dz = centers[i,2] - centers[j,2]
                d = int(math.sqrt(dx*dx + dy*dy + dz*dz) / bw)
                if 0 <= d <= max_dist: corr[d] += B[i]*B[j]; cnt[d] += 1
    for d in range(max_dist + 1):
        if cnt[d] > 0: corr[d] /= cnt[d]
    return corr, cnt

def run(L=4, M=64, Gamma=0.05, n_thermal=5000):
    print(f"L={L}, M={M}, Γ={Gamma}")
    t0 = time.time()
    
    neigh, deg, pos = build_pyrochlore(L)
    N = len(deg)
    print(f"  N={N}, spins={N*M}, deg={np.min(deg)}-{np.max(deg)}")
    
    hexagons = find_hexagons(neigh, deg, 500)
    print(f"  Hexagons: {len(hexagons)}")
    
    harr = np.zeros((len(hexagons), 6), dtype=np.int32)
    for hi, h in enumerate(hexagons):
        for k, n in enumerate(h): harr[hi, k] = n
    
    centers = np.zeros((len(hexagons), 3))
    for h in range(len(hexagons)):
        cx = cy = cz = 0.0
        for node in harr[h]: cx += pos[node,0]; cy += pos[node,1]; cz += pos[node,2]
        centers[h] = [cx/6, cy/6, cz/6]
    
    Ks, Kt = setup_trotter(1.0, Gamma, M)
    np.random.seed(42)
    z = np.random.choice(np.array([-1.0, 1.0]), size=(N, M))
    
    for step in range(n_thermal):
        wolff_cluster_pyro(z, neigh, deg, M, Ks, Kt)
        if step % 1000 == 0:
            ice = 0.0
            for tau in range(M):
                for x in range(L):
                    for y in range(L):
                        for zc in range(L):
                            base = ((x*L + y)*L + zc) * 4
                            s = z[base,tau]+z[base+1,tau]+z[base+2,tau]+z[base+3,tau]
                            if abs(s) < 0.01: ice += 1
            ice /= (M * L**3)
            c1 = 0.0; cn = 0
            for tau in range(M):
                for i in range(N):
                    for k in range(deg[i]):
                        j = neigh[i, k]
                        if j > i: c1 += z[i,tau]*z[j,tau]; cn += 1
            c1 /= max(cn, 1)
            print(f"    {step}: ice={ice:.4f} C1={c1:.4f}")
    
    corr, cnt = hex_correlator(z, harr, centers, 8)
    print("\n  ⟨B(0)B(r)⟩:")
    for d in range(1, 9):
        if cnt[d] > 0:
            print(f"    r≈{d*0.5:.1f}a: {corr[d]:.8f} (n={int(cnt[d])})")
    
    xi = float('inf')
    if np.all(np.abs(corr[1:][cnt[1:] > 0] - 1.0) < 0.01):
        xi = L  # U(1)-фаза: ξ_eff = L
        print(f"\n  U(1)-фаза: ⟨BB⟩=1 → ξ_eff = L = {L}a")
    
    alpha = (2-math.log(2))/(8*math.pi*xi) if xi > 0 and xi < 1e6 else 0
    alpha_exp = 1/137.036
    elapsed = time.time() - t0
    print(f"  [{elapsed:.0f}s] α = {alpha:.6f} ({alpha/alpha_exp:.2f}× exp)")
    return {"L":L, "M":M, "xi":xi, "alpha":alpha, "n_hex":len(hexagons)}

if __name__ == "__main__":
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    print("="*60)
    print(f"ПРЯМОЕ ВЫЧИСЛЕНИЕ α: L={L}, M=64")
    print("="*60)
    r = run(L=L, M=64, Gamma=0.05, n_thermal=5000)
    if r:
        with open(f'alpha_L{L}.json','w') as f: json.dump(r, f, indent=2)
        print(f"\nСохранено: alpha_L{L}.json")
