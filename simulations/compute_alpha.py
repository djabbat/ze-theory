"""
Вычисление α из принципов Ze.
Формула: α = P(T|v*)/(4π·ξ/a) = (2-ln2)/(8π·ξ/a)
"""
import math
v_star = 1 - math.log(2)
PT = (2 - math.log(2)) / 2
alpha_exp = 1/137.036

print("=" * 55)
print("ВЫЧИСЛЕНИЕ α ИЗ ПРИНЦИПОВ Ze")
print("=" * 55)
print(f"v* = 1 - ln 2 = {v_star:.6f}")
print(f"P(T|v*) = {PT:.4f}")
print(f"Формула: α = P(T|v*)/(4π · ξ/a) = {PT:.4f}/(4π · ξ/a)")
print()

for xi in [3,4,5,6,7,8,10,16,32]:
    alpha = PT/(4*math.pi*xi) if xi > 0 else 0
    ratio = alpha/alpha_exp
    print(f"  ξ={xi:2}a: α={alpha:.6f} ({ratio:.2f}× exp)")

xi_exact = PT/(4*math.pi*alpha_exp)
print(f"\n  ξ={xi_exact:.2f}a: α={alpha_exp:.8f} ← ТОЧНО 1/137")
print(f"\nДля экспериментального α: ξ/a = {xi_exact:.2f}")
print(f"Фундаментальный масштаб Ze: a = λ_e/{xi_exact:.2f} ≈ {3.86e-13/xi_exact:.1e} м")
