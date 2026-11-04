# Для проверки законов сохранения
from math import pi


def p(m1, m2, v1, v2):
    return m1 * v1 + m2 * v2


def Ek(m1, m2, vx1, vy1, vx2, vy2):
    return (m1 * (vx1**2 + vy1**2) + m2 * (vx2**2 + vy2**2)) / 2


def p_ort(m, vx, vy, nx, ny):
    return m * (-vx * ny + vy * nx)


# Максимальная деформация и длительность лобового контакта
# u - относительная скорость сближения, mu = m1*m2/(m1 + m2)


def delta_max(mu, u, k, mode):
    if mode == 1:
        return u * (mu / k) ** 0.5
    return (5 * mu * u**2 / (4 * k)) ** 0.4


def t_contact(mu, u, k, mode):
    if mode == 1:
        return pi * (mu / k) ** 0.5
    return 2.9432752 * delta_max(mu, u, k, mode) / u
