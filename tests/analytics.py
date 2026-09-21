# Для проверки законов сохранения


def p(m1, m2, v1, v2):
    return m1 * v1 + m2 * v2


def Ek(m1, m2, vx1, vy1, vx2, vy2):
    return (m1 * (vx1**2 + vy1**2) + m2 * (vx2**2 + vy2**2)) / 2


def p_ort(m, vx, vy, nx, ny):
    return m * (-vx * ny + vy * nx)
