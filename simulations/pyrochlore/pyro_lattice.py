"""Общий модуль: геометрия пирохлорной решётки (координация 6)."""
import numpy as np
from numba import njit

@njit
def build_pyrochlore(L):
    """
    Правильная пирохлорная решётка. N = 4L³ спинов.
    up-тетраэдр + down-тетраэдр для каждой примитивной ячейки.
    Каждый спин имеет ровно 6 соседей.
    
    Возвращает: neigh[N,6], deg[N], pos[N,3]
    """
    N = 4 * L**3
    max_n = 6
    neigh = -np.ones((N, max_n), dtype=np.int32)
    deg = np.zeros(N, dtype=np.int32)
    pos = np.zeros((N, 3), dtype=np.float64)
    
    def idx(x, y, z, s):
        return (((x % L) * L + (y % L)) * L + (z % L)) * 4 + s
    
    def add(i, j):
        if 0 <= i < N and 0 <= j < N and deg[i] < max_n:
            for k in range(deg[i]):
                if neigh[i, k] == j:
                    return
            neigh[i, deg[i]] = j
            deg[i] += 1
    
    # Базисные векторы подрешёток (в единицах кубической ячейки)
    basis = [(0, 0, 0), (0.25, 0.25, 0), (0.25, 0, 0.25), (0, 0.25, 0.25)]
    
    for x in range(L):
        for y in range(L):
            for z in range(L):
                # Позиции 4 спинов в ячейке
                for s, (dx, dy, dz) in enumerate(basis):
                    pos[idx(x, y, z, s)] = [x + dx, y + dy, z + dz]
                
                # UP-тетраэдр: все 4 спина в одной ячейке
                su = [idx(x, y, z, s) for s in range(4)]
                for a in range(4):
                    for b in range(a + 1, 4):
                        add(su[a], su[b])
                        add(su[b], su[a])
                
                # DOWN-тетраэдр: s0(x,y,z), s1(x-1,y-1,z), s2(x-1,y,z-1), s3(x,y-1,z-1)
                sd = [
                    idx(x, y, z, 0),
                    idx(x - 1, y - 1, z, 1),
                    idx(x - 1, y, z - 1, 2),
                    idx(x, y - 1, z - 1, 3),
                ]
                for a in range(4):
                    for b in range(a + 1, 4):
                        add(sd[a], sd[b])
                        add(sd[b], sd[a])
    
    return neigh, deg, pos


def test_lattice(L=3):
    """Проверка корректности построения решётки."""
    neigh, deg, pos = build_pyrochlore(L)
    N = len(deg)
    dmin, dmax = np.min(deg), np.max(deg)
    
    assert dmin == 6, f"Минимальная степень {dmin}, ожидается 6"
    assert dmax == 6, f"Максимальная степень {dmax}, ожидается 6"
    assert N == 4 * L**3, f"N={N}, ожидается {4*L**3}"
    
    # Проверка симметричности
    for i in range(N):
        for k in range(deg[i]):
            j = neigh[i, k]
            found = False
            for k2 in range(deg[j]):
                if neigh[j, k2] == i:
                    found = True
                    break
            assert found, f"Связь {i}→{j} не симметрична"
    
    return True


if __name__ == "__main__":
    for L in [2, 3, 4, 5]:
        ok = test_lattice(L)
        print(f"L={L}: {'✅' if ok else '❌'}")
    print("Все тесты пройдены!")
