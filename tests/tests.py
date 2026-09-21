from math import cos, sin, radians
from src.main import model_contact
from tests.analytics import momentum, kinetic_energy, transverse_momentum


total = 0
passed = 0


def check(name, got, expected, max_diff):
    global total, passed
    total += 1
    diff = abs(got - expected)
    if diff <= max_diff:
        passed += 1
    status = "OK" if diff <= max_diff else "WRONG"
    print(f"  [{status}] {name:<24} {got:>12.6f} {expected:>12.6f} {diff:>10.2e}")


def test(mode, m1, m2, u1, u2, angle=0):
    nx, ny = cos(radians(angle)), sin(radians(angle))
    ux1, uy1 = u1[0]*nx - u1[1]*ny, u1[0]*ny + u1[1]*nx
    ux2, uy2 = u2[0]*nx - u2[1]*ny, u2[0]*ny + u2[1]*nx
    vx1, vx2, vy1, vy2, *_ = model_contact(
        1, 1, ux1, ux2, uy1, uy2, 0, 2*nx, 0, 2*ny, m1, m2, 1e8, mode)

    print(f"\n{'Гук' if mode == 1 else 'Герц'}: m = ({m1}, {m2}), u = {u1}, {u2}, угол = {angle}°")
    check("импульс x", momentum(m1, m2, vx1, vx2), momentum(m1, m2, ux1, ux2), 1e-8)
    check("импульс y", momentum(m1, m2, vy1, vy2), momentum(m1, m2, uy1, uy2), 1e-8)

    e0 = kinetic_energy(m1, m2, ux1, uy1, ux2, uy2)
    check("кинетическая энергия", kinetic_energy(m1, m2, vx1, vy1, vx2, vy2), e0, 1e-7*e0)

    check("поперечный импульс 1", transverse_momentum(m1, vx1, vy1, nx, ny),
          transverse_momentum(m1, ux1, uy1, nx, ny), 1e-3)
    check("поперечный импульс 2", transverse_momentum(m2, vx2, vy2, nx, ny),
          transverse_momentum(m2, ux2, uy2, nx, ny), 1e-3)


if __name__ == "__main__":
    for mode in (1, 2):
        # Лобовые удары: равные и разные массы, неподвижные и движущиеся шары.
        test(mode, 1, 1, (2, 0), (0, 0))
        test(mode, 2, 1, (2, 0), (0, 0))
        test(mode, 1, 2, (2, 0), (0, 0))
        test(mode, 1, 1, (2, 0), (-1, 0))
        test(mode, 2, 1, (2, 0), (-1, 0))
        test(mode, 1, 1, (3, 0), (1, 0))  # Первый догоняет второй
        test(mode, 1, 1, (0, 0), (-2, 0))  # Первый шар покоится

        # Косые удары
        test(mode, 1, 1, (2, 1), (0, 0))
        test(mode, 2, 1, (2, 1), (-1, 0.5))
        test(mode, 1, 1, (2, -1), (0, 0))
        test(mode, 2, 1, (2, 1), (-1, 0.5), 90)
        test(mode, 1, 1, (2, 1), (0, 0), 45)

    print(f"\n{passed} / {total}")
    raise SystemExit(0 if passed == total else 1)
