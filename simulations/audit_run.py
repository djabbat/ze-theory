"""
Аудит-версия симуляторов Ze. Все баги исправлены.
Запуск: python3 audit_run.py
"""
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

# ============================================================
# 1. Классическое MC с staggered magnetization
# ============================================================
def run_classical_scan():
    from classical_mc.ze_mc import run_metropolis, total_energy
    np.random.seed(42)
    print("=" * 70)
    print("КЛАССИЧЕСКОЕ MC — фазовая диаграмма J_s vs T")
    print(f"{'T':>8} {'J_s':>8} {'|v|':>10} {'|v_stag|':>10} {'E/N':>10}")
    print("-" * 50)
    
    Lx, Ly, Lt = 4, 4, 8
    J_t, h = 1.0, 0.0
    
    for T in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
        for J_s in [0.0, 0.1, 0.2, 0.3, 0.5]:
            beta = 1.0 / T
            z = np.random.choice([-1, 1], size=(Lx, Ly, Lt)).astype(np.float64)
            
            _, _, _ = run_metropolis(z, Lx, Ly, Lt, J_t, J_s, h, beta, 2000, 50)
            
            v = np.mean(z)
            even = np.mean(z[:,:,0::2]); odd = np.mean(z[:,:,1::2])
            v_stag = abs(even - odd)
            E = total_energy(z, Lx, Ly, Lt, J_t, J_s, h) / (Lx*Ly*Lt)
            
            phase = "АФМ" if v_stag > 0.3 else ("ФМ" if abs(v) > 0.3 else "пара")
            print(f"{T:8.2f} {J_s:8.2f} {abs(v):10.4f} {v_stag:10.4f} {E:10.4f}  {phase}")
    print()

# ============================================================
# 2. Квантовое MC с правильным знаком
# ============================================================
def run_quantum_scan():
    print("=" * 70)
    print("КВАНТОВОЕ MC — сканирование по Γ")
    print(f"{'Γ':>8} {'|v|':>10} {'|v_stag|':>10} {'E/N':>10}")
    print("-" * 50)
    
    from quantum_mc.ze_qmc import QuantumZeSimulation
    
    for Gamma in [0.2, 0.5, 0.8, 1.0, 1.5, 2.0]:
        sim = QuantumZeSimulation(L=8, M_trotter=16, Gamma=Gamma, J_t=1.0,
                                   n_thermal=1000, n_samples=3000)
        r = sim.run()
        z = sim.final_config
        even = np.mean(z[0::2,:]); odd = np.mean(z[1::2,:])
        v_stag = abs(even - odd)
        phase = "АФМ" if v_stag > 0.2 else "пара"
        print(f"{Gamma:8.2f} {r['mean_v']:10.4f} {v_stag:10.4f} {r['mean_energy']:10.4f}  {phase}")
    print()

# ============================================================
if __name__ == "__main__":
    run_classical_scan()
    run_quantum_scan()
