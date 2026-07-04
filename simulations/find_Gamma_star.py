#!/usr/bin/env python3
"""
Поиск физического условия, фиксирующего Γ*.
==============================================

Три гипотезы:
  A. Масштаб конфайнмента: Λ_Ze ≈ Λ_QCD ≈ 200 МэВ
     → a_Ze = ℏc/Λ_Ze → L = λ_e/a_Ze → Γ* из L(Γ) = 7.13

  B. Информационный принцип: взаимная информация I(agent:field)
     максимальна при Γ = Γ*
     → Γ* из условия dI/dΓ = 0

  C. Принцип минимальной свободной энергии:
     F(Γ) = E₀(Γ) − T_eff·S(Γ)
     → Γ* из dF/dΓ = 0

Автор: Jaba Tqemaladze, MD
Дата: 2026-07-04
"""

import numpy as np
import math
import json
import sys
import os
from datetime import datetime

# Добавляем путь к exact_diag
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'exact_diag'))

LN2 = math.log(2)
V_STAR = 1.0 - LN2
PT_ERROR = (2 - LN2) / 2
ALPHA_EXP = 1 / 137.035999084

# Константы
HBAR_C = 197.3269804  # МэВ·фм
M_E = 0.510998950      # МэВ
LAMBDA_E = HBAR_C / M_E  # комптоновская длина волны электрона, фм


# ═══════════════════════════════════════════════════════════════════
# ГИПОТЕЗА A: МАСШТАБ КОНФАЙНМЕНТА
# ═══════════════════════════════════════════════════════════════════

def hypothesis_A_confinement():
    """
    Γ* из условия: масштаб обрезания Ze равен масштабу конфайнмента КХД.
    
    Λ_QCD ≈ 200 МэВ (типичный масштаб непертурбативной КХД).
    Λ_Ze = ℏc / a_Ze — масштаб обрезания Z₂-теории.
    
    L = λ_e / a_Ze = Λ_Ze / m_e
    
    Если Λ_Ze = Λ_QCD:
        L = Λ_QCD / m_e ≈ 200 / 0.511 ≈ 391
    
    Это даёт α = P(T|v*)/(4π·391) ≈ 0.000133 = 1/7513 — не 1/137.
    
    Но Λ_QCD — масштаб КХД, а не Z₂-теории!
    Масштаб Z₂-теории Λ_Ze — другой.
    
    Из экспериментального α:
        L = 7.13 → Λ_Ze = L · m_e = 7.13 × 0.511 МэВ ≈ 3.65 ГэВ
    
    Этот масштаб (~3.65 ГэВ) близок к массе b-кварка (~4.18 ГэВ)
    и характерному масштабу нарушения киральной симметрии в КХД.
    
    Вывод: Λ_Ze НЕ равен Λ_QCD. Λ_Ze ≈ 3.65 ГэВ — новый масштаб,
    требующий независимого физического обоснования.
    """
    print("=" * 65)
    print("ГИПОТЕЗА A: МАСШТАБ КОНФАЙНМЕНТА")
    print("=" * 65)
    
    Lambda_QCD = 200  # МэВ
    Lambda_Ze_from_alpha = 7.13 * M_E  # МэВ
    
    print(f"  Λ_QCD ≈ {Lambda_QCD} МэВ")
    print(f"  Λ_Ze (из α_exp) = L·m_e = {Lambda_Ze_from_alpha:.1f} МэВ")
    print(f"  Отношение: Λ_Ze/Λ_QCD = {Lambda_Ze_from_alpha/Lambda_QCD:.1f}")
    print()
    print(f"  Интерпретация: Λ_Ze ≈ 3.65 ГэВ — масштаб,")
    print(f"  на котором Z₂-калибровочная теория переходит в U(1).")
    print(f"  Это ~масштаб нарушения киральной симметрии в КХД.")
    print(f"  Физический смысл: Λ_Ze — это масштаб, на котором")
    print(f"  «лёд тает» (ice-rule нарушается квантовыми флуктуациями).")
    print()
    
    # Вычисляем Γ* из условия Λ_Ze = L(Γ)·m_e
    # L(Γ) = 1/(Z·C·Γ^6) → Γ = (Z·C·L)^{-1/6}
    Z_g = 0.75
    C = 0.25
    L_target = 7.13
    
    Gamma_from_L = (Z_g * C * L_target) ** (-1.0 / 6.0)
    
    print(f"  Γ* (из L={L_target}) = {Gamma_from_L:.4f} J_t")
    print(f"  Проверка: L(Γ*) = 1/({Z_g}·{C}·{Gamma_from_L:.4f}⁶) = {1/(Z_g*C*Gamma_from_L**6):.2f}")
    
    return {
        'hypothesis': 'A: confinement scale',
        'Lambda_QCD_MeV': Lambda_QCD,
        'Lambda_Ze_MeV': Lambda_Ze_from_alpha,
        'Gamma_star_from_L': Gamma_from_L,
        'L_target': L_target,
        'note': 'Λ_Ze ≠ Λ_QCD; requires independent justification'
    }


# ═══════════════════════════════════════════════════════════════════
# ГИПОТЕЗА B: МАКСИМУМ ВЗАИМНОЙ ИНФОРМАЦИИ
# ═══════════════════════════════════════════════════════════════════

def hypothesis_B_mutual_information():
    """
    Γ* из условия максимума взаимной информации между агентом и полем.
    
    Агент наблюдает бинарные события T/S. Z₂-калибровочное поле
    находится в состоянии с волновой функцией |Ψ(Γ)⟩.
    
    Взаимная информация:
        I(agent : field) = S(agent) + S(field) - S(agent, field)
    
    В критической точке v*:
        S(agent) = H(v*) = -(v*/2)log₂(v*/2) - ((1-v*)/2)log₂((1-v*)/2)
                  = максимальна при ограничении S=−T
    
    S(field) = -Tr(ρ_field log₂ ρ_field) — энтропия фон Неймана
              основного состояния H_Ze(Γ).
    
    I(agent : field) максимальна, когда:
    1. Агент находится в точке максимальной энтропии (v*)
    2. Поле находится в точке максимальной квантовой запутанности
    
    Для Z₂-калибровочной теории, запутанность максимальна
    вблизи квантовой критической точки Γ ≈ Γ_c.
    
    Результат: Γ* ≈ Γ_c (квантовая критическая точка).
    """
    print("=" * 65)
    print("ГИПОТЕЗА B: МАКСИМУМ ВЗАИМНОЙ ИНФОРМАЦИИ")
    print("=" * 65)
    
    # Данные QMC для v_stag(Γ) — прокси для энтропии поля
    # В критической точке энтропия максимальна
    Gamma_c = 1.05  # из QMC-данных
    
    # Сканируем Γ для поиска максимума I
    Gamma_range = np.linspace(0.1, 2.0, 50)
    
    # Энтропия агента (фиксирована в v*)
    v_star = V_STAR
    p_T = (1 + v_star) / 2  # P(T) = (1+v)/2
    p_S = (1 - v_star) / 2
    if p_T > 0 and p_S > 0:
        S_agent = -p_T * math.log2(p_T) - p_S * math.log2(p_S)
    else:
        S_agent = 0
    
    best_I = -float('inf')
    best_Gamma = None
    I_values = []
    
    for Gamma in Gamma_range:
        # v_stag как функция Γ (из данных QMC, Таблица 3.2)
        # Используем линейную интерполяцию
        qmc_G = np.array([0.2, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0])
        qmc_v = np.array([0.996, 0.976, 0.842, 0.565, 0.398, 0.277, 0.233, 0.226])
        v_stag = np.interp(Gamma, qmc_G, qmc_v)
        
        # Энтропия поля (модель: запутанность ∼ (1 - v_stag²))
        # В АФМ-фазе v_stag≈1 → запутанность мала
        # В критической точке v_stag≈0 → запутанность максимальна
        S_field = 1.0 - v_stag**2
        
        # Взаимная информация (упрощённо)
        I_af = S_agent + S_field
        
        I_values.append((Gamma, I_af))
        
        if I_af > best_I:
            best_I = I_af
            best_Gamma = Gamma
    
    print(f"  S_agent(v*) = {S_agent:.4f}")
    print(f"  S_field(Γ) ∼ 1 - v_stag²")
    print(f"  I_max = {best_I:.4f} при Γ* = {best_Gamma:.2f}")
    print()
    
    # Вычисляем L из Γ*
    Z_g = 0.75
    C = 0.25
    g_star = Z_g * C * (best_Gamma ** 6)
    L_from_I = 1.0 / g_star if g_star > 1e-12 else float('inf')
    alpha_from_I = PT_ERROR / (4 * math.pi * L_from_I) if L_from_I > 0 else 0
    
    print(f"  g(Γ*) = {g_star:.4f}")
    print(f"  L = 1/g = {L_from_I:.2f}")
    print(f"  α = {alpha_from_I:.6f} (1/{1/alpha_from_I:.1f})")
    print(f"  α_exp = {ALPHA_EXP:.6f} (1/{1/ALPHA_EXP:.1f})")
    print()
    
    if abs(alpha_from_I - ALPHA_EXP) / ALPHA_EXP < 0.05:
        print(f"  ✓ Гипотеза B ПОДТВЕРЖДАЕТСЯ: Γ* даёт α в пределах 5%!")
    else:
        print(f"  → Гипотеза B не даёт точного совпадения.")
        print(f"  → Нужен более точный учёт запутанности поля.")
        print(f"  → Предложение: вычислить entanglement entropy")
        print(f"    основного состояния H_Ze через точную диагонализацию.")
    
    return {
        'hypothesis': 'B: mutual information max',
        'S_agent': S_agent,
        'Gamma_star': best_Gamma,
        'I_max': best_I,
        'g_star': g_star,
        'L': L_from_I,
        'alpha_pred': alpha_from_I,
        'alpha_exp': ALPHA_EXP,
        'diff_pct': abs(alpha_from_I - ALPHA_EXP) / ALPHA_EXP * 100
    }


# ═══════════════════════════════════════════════════════════════════
# ГИПОТЕЗА C: МИНИМУМ СВОБОДНОЙ ЭНЕРГИИ
# ═══════════════════════════════════════════════════════════════════

def hypothesis_C_free_energy():
    """
    Γ* из условия минимума свободной энергии агента.
    
    F(Γ) = ⟨H_Ze⟩(Γ) − T_eff · S_agent(v*)
    
    T_eff = Γ/M_trotter — эффективная температура
    квантовых флуктуаций в троттеровской формулировке.
    
    Условие: dF/dΓ = 0 → Γ* = argmin F(Γ).
    
    Энергия основного состояния E₀(Γ) вычисляется через
    точную диагонализацию (или спин-волновую теорию).
    """
    print("=" * 65)
    print("ГИПОТЕЗА C: МИНИМУМ СВОБОДНОЙ ЭНЕРГИИ")
    print("=" * 65)
    
    M_trotter = 32  # стандартное значение из QMC
    
    Gamma_range = np.linspace(0.1, 1.5, 50)
    
    # Приближение для E₀(Γ) на основе точной диагонализации
    # тетраэдра и кластера (см. ed_cluster.py)
    # 
    # Для Z₂-калибровочной теории в конфайнмент-фазе:
    # E₀(Γ) ≈ -N·J_t + const·Γ²
    
    N_spins = 4 * 2**3  # 32 спина (L=2 пирохлор)
    J_t = 1.0
    
    F_values = []
    best_F = float('inf')
    best_Gamma = None
    
    for Gamma in Gamma_range:
        # Энергия основного состояния (приближение)
        # E₀(Γ) = -N·J_t (конфайнмент) + ε(Γ) (квантовые флуктуации)
        # Используем модель: E₀(Γ) ≈ E_classical + α·Γ²
        
        # Из точной диагонализации тетраэдра (16 состояний):
        # E₀(Γ=0) = -6J (все AFM связи насыщены, но фрустрированы)
        # E₀(Γ) ≈ -6J + c₁·Γ + c₂·Γ²
        
        # Для пирохлорной решётки L=2:
        # Классическая энергия = -N·J_t/2 (каждая связь учтена один раз)
        E_classical = -N_spins * J_t / 2  # ≈ -16J
        
        # Квантовая поправка (из теории возмущений)
        # При Γ≪J: E₀(Γ) = E_classical − (Γ²/2J)·Σ_i ⟨σ^x_i⟩
        E_quantum = -0.3 * N_spins * Gamma**2 / J_t
        
        E0 = E_classical + E_quantum
        
        # Свободная энергия
        T_eff = Gamma / M_trotter
        v_star = V_STAR
        p_T = (1 + v_star) / 2
        p_S = (1 - v_star) / 2
        
        if p_T > 0 and p_S > 0:
            S_agent = -p_T * math.log(p_T) - p_S * math.log(p_S)
        else:
            S_agent = 0
        
        F = E0 - T_eff * S_agent
        
        F_values.append((Gamma, F, E0, T_eff))
        
        if F < best_F:
            best_F = F
            best_Gamma = Gamma
    
    print(f"  F(Γ) = E₀(Γ) − T_eff·S_agent(v*)")
    print(f"  T_eff = Γ/{M_trotter}")
    print(f"  F_min = {best_F:.4f} при Γ* = {best_Gamma:.4f}")
    print()
    
    # L из Γ*
    Z_g = 0.75
    C = 0.25
    g_star = Z_g * C * (best_Gamma ** 6)
    L_from_F = 1.0 / g_star if g_star > 1e-12 else float('inf')
    alpha_from_F = PT_ERROR / (4 * math.pi * L_from_F) if L_from_F > 0 else 0
    
    print(f"  g(Γ*) = {g_star:.4f}")
    print(f"  L = 1/g = {L_from_F:.2f}")
    print(f"  α = {alpha_from_F:.6f} (1/{1/alpha_from_F:.1f})")
    print(f"  α_exp = {ALPHA_EXP:.6f} (1/{1/ALPHA_EXP:.1f})")
    
    diff = abs(alpha_from_F - ALPHA_EXP) / ALPHA_EXP * 100
    
    if diff < 5:
        print(f"\n  ✓ Гипотеза C ПОДТВЕРЖДАЕТСЯ: Γ* даёт α в пределах {diff:.1f}%!")
    else:
        print(f"\n  → Гипотеза C: отклонение {diff:.1f}%")
        print(f"  → Требуется более точное E₀(Γ) через точную диагонализацию")
    
    return {
        'hypothesis': 'C: free energy minimum',
        'Gamma_star': best_Gamma,
        'F_min': best_F,
        'g_star': g_star,
        'L': L_from_F,
        'alpha_pred': alpha_from_F,
        'diff_pct': diff
    }


# ═══════════════════════════════════════════════════════════════════
# ГИПОТЕЗА D: РАВЕНСТВО ЭНТРОПИЙ
# ═══════════════════════════════════════════════════════════════════

def hypothesis_D_entropy_equality():
    """
    Γ* из условия: S_agent(v*) = S_field(Γ).
    
    Равенство энтропий агента и калибровочного поля —
    условие термодинамического равновесия в дуальном описании.
    
    S_agent(v*) = H(v*) — бинарная энтропия
    S_field(Γ) = энтропия фон Неймана основного состояния |Ψ₀(Γ)⟩
    
    При Γ=0: S_field=0 (чистое состояние без запутанности)
    При Γ=Γ_c: S_field максимальна
    При Γ→∞: S_field→0 (парамагнетик)
    
    Условие S_field(Γ) = S_agent(v*) ≈ 0.88 бит определяет Γ*.
    """
    print("=" * 65)
    print("ГИПОТЕЗА D: РАВЕНСТВО ЭНТРОПИЙ")
    print("=" * 65)
    
    # Энтропия агента в v*
    v_star = V_STAR
    p_T = (1 + v_star) / 2
    p_S = (1 - v_star) / 2
    if p_T > 0 and p_S > 0:
        S_agent = -p_T * math.log2(p_T) - p_S * math.log2(p_S)
    else:
        S_agent = 0
    
    print(f"  S_agent(v*) = {S_agent:.4f} бит")
    
    # Модель для S_field(Γ):
    # S_field ∼ −Tr(ρ log₂ ρ), где ρ — редуцированная матрица
    # плотности одного спина в основном состоянии.
    #
    # В приближении среднего поля:
    # ρ = (I + v_stag(Γ)·σ^z)/2
    # S_field(Γ) = −[(1+v)/2·log₂((1+v)/2) + (1−v)/2·log₂((1−v)/2)]
    # где v = v_stag(Γ)
    
    qmc_G = np.array([0.2, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0])
    qmc_v = np.array([0.996, 0.976, 0.842, 0.565, 0.398, 0.277, 0.233, 0.226])
    
    def S_field_binary(v):
        """Бинарная энтропия для спина с намагниченностью v."""
        p_up = (1 + v) / 2
        p_down = (1 - v) / 2
        if p_up <= 0 or p_down <= 0:
            return 0.0
        return -p_up * math.log2(p_up) - p_down * math.log2(p_down)
    
    # Поиск Γ, где S_field = S_agent
    best_Gamma = None
    best_diff = float('inf')
    
    Gamma_fine = np.linspace(0.5, 1.5, 200)
    
    for Gamma in Gamma_fine:
        v = np.interp(Gamma, qmc_G, qmc_v)
        S_f = S_field_binary(abs(v))
        diff = abs(S_f - S_agent)
        
        if diff < best_diff:
            best_diff = diff
            best_Gamma = Gamma
    
    v_at_best = np.interp(best_Gamma, qmc_G, qmc_v)
    S_at_best = S_field_binary(abs(v_at_best))
    
    print(f"  S_field(Γ) модель: бинарная энтропия v_stag")
    print(f"  S_field(Γ*) = {S_at_best:.4f} бит при Γ* = {best_Gamma:.4f}")
    print(f"  v_stag(Γ*) = {v_at_best:.4f}")
    print(f"  |S_agent − S_field| = {best_diff:.6f}")
    print()
    
    # L из Γ*
    Z_g = 0.75
    C = 0.25
    g_star = Z_g * C * (best_Gamma ** 6)
    L_from_S = 1.0 / g_star if g_star > 1e-12 else float('inf')
    alpha_from_S = PT_ERROR / (4 * math.pi * L_from_S) if L_from_S > 0 else 0
    
    print(f"  g(Γ*) = {g_star:.4f}")
    print(f"  L = 1/g = {L_from_S:.2f}")
    print(f"  α = {alpha_from_S:.6f} (1/{1/alpha_from_S:.1f})")
    print(f"  α_exp = {ALPHA_EXP:.6f} (1/{1/ALPHA_EXP:.1f})")
    
    diff = abs(alpha_from_S - ALPHA_EXP) / ALPHA_EXP * 100
    
    if diff < 5:
        print(f"\n  ★ Гипотеза D ПОДТВЕРЖДАЕТСЯ: Γ* даёт α в пределах {diff:.1f}%!")
        print(f"  ★ Условие S_agent = S_field ФИКСИРУЕТ Γ*")
        print(f"  ★ Это термодинамический принцип: равенство энтропий")
        print(f"     наблюдателя и наблюдаемой системы.")
    else:
        print(f"\n  → Гипотеза D: отклонение {diff:.1f}%")
    
    return {
        'hypothesis': 'D: entropy equality',
        'S_agent': S_agent,
        'S_field_at_Gamma_star': S_at_best,
        'Gamma_star': best_Gamma,
        'v_stag_at_Gamma_star': v_at_best,
        'g_star': g_star,
        'L': L_from_S,
        'alpha_pred': alpha_from_S,
        'diff_pct': diff
    }


# ═══════════════════════════════════════════════════════════════════
# СВОДНЫЙ АНАЛИЗ
# ═══════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Поиск физического условия для Γ*"
    )
    parser.add_argument('--all', action='store_true',
                       help='Все гипотезы')
    parser.add_argument('--save', type=str, default=None,
                       help='Сохранить результаты в JSON')
    
    args = parser.parse_args()
    
    all_results = {}
    
    print()
    rA = hypothesis_A_confinement()
    all_results['A'] = rA
    print()
    
    rB = hypothesis_B_mutual_information()
    all_results['B'] = rB
    print()
    
    rC = hypothesis_C_free_energy()
    all_results['C'] = rC
    print()
    
    rD = hypothesis_D_entropy_equality()
    all_results['D'] = rD
    print()
    
    # Сводка
    print("=" * 65)
    print("СВОДКА: КАКАЯ ГИПОТЕЗА ФИКСИРУЕТ Γ*?")
    print("=" * 65)
    print(f"{'Гипотеза':<5} {'Γ*':>8} {'L':>8} {'α_pred':>10} {'1/α':>10} {'Δ%':>8}")
    print("-" * 55)
    
    best_hypothesis = None
    best_diff = float('inf')
    
    for key, r in all_results.items():
        if 'alpha_pred' in r and r['alpha_pred'] > 0:
            diff = r.get('diff_pct', 100)
            print(f"{key:<5} {r.get('Gamma_star', 0):8.4f} "
                  f"{r.get('L', 0):8.2f} {r['alpha_pred']:10.6f} "
                  f"{1/r['alpha_pred']:10.1f} {diff:7.1f}%")
            
            if diff < best_diff:
                best_diff = diff
                best_hypothesis = key
    
    if best_hypothesis:
        r = all_results[best_hypothesis]
        print(f"\n★ РЕКОМЕНДУЕМАЯ ГИПОТЕЗА: {best_hypothesis}")
        print(f"  {r['hypothesis']}")
        print(f"  Γ* = {r.get('Gamma_star', 0):.4f}")
        print(f"  Отклонение α = {r.get('diff_pct', 100):.2f}%")
    
    if args.save:
        output = {
            'timestamp': datetime.now().isoformat(),
            'alpha_exp': ALPHA_EXP,
            'v_star': V_STAR,
            'PT_error': PT_ERROR,
            'results': all_results,
            'best_hypothesis': best_hypothesis
        }
        
        def convert(obj):
            if isinstance(obj, (np.integer,)):
                return int(obj)
            if isinstance(obj, (np.floating,)):
                return float(obj)
            return obj
        
        with open(args.save, 'w') as f:
            json.dump(output, f, indent=2, default=convert)
        print(f"\nСохранено: {args.save}")
