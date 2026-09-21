"""Таблицы отчёта: python docs/experiments.py (окно графика не открывается)."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
# Расчёт идеально упругого мгновенного удара
from math import hypot


def elastic_velocities(m1, m2, u1, u2, p1, p2):
    # Сохраняем касательные скорости, обращаем нормальную относительную
    dx, dy = p1[0] - p2[0], p1[1] - p2[1]
    length = hypot(dx, dy)
    nx, ny = dx / length, dy / length
    un = (u1[0] - u2[0]) * nx + (u1[1] - u2[1]) * ny
    impulse = -2 * un / (1 / m1 + 1 / m2)
    return ((u1[0] + impulse * nx / m1, u1[1] + impulse * ny / m1),
            (u2[0] - impulse * nx / m2, u2[1] - impulse * ny / m2))


def momentum(m1, m2, v1, v2):
    return tuple(m1 * a + m2 * b for a, b in zip(v1, v2))


def energy(m1, m2, v1, v2):
    return (m1 * sum(v * v for v in v1) + m2 * sum(v * v for v in v2)) / 2


# Начальные условия в момент первого касания
from src.main import model_contact

# Название, массы, скорости, координаты центров. Радиусы обоих шаров — 1 м.
CASES = [
    ("Лобовой, равные массы, второй покоится", 1., 1., (2., 0.), (0., 0.), (0., 0.), (2., 0.)),
    ("Лобовой, разные массы, второй покоится", 2., 1., (2., 0.), (0., 0.), (0., 0.), (2., 0.)),
    ("Лобовой, равные массы, встречное движение", 1., 1., (2., 0.), (-1., 0.), (0., 0.), (2., 0.)),
    ("Лобовой, разные массы, встречное движение", 2., 1., (2., 0.), (-1., 0.), (0., 0.), (2., 0.)),
    ("Косой, равные массы, второй покоится", 1., 1., (2., 1.), (0., 0.), (0., 0.), (2., 0.)),
    ("Косой, разные массы, оба движутся", 2., 1., (2., 1.), (-1., .5), (0., 0.), (2., 0.)),
    ("Лобовой вдоль оси y", 1., 1., (0., 2.), (0., 0.), (0., 0.), (0., 2.)),
]


def simulate_case(case, mode, k=1e8):
    _, m1, m2, u1, u2, p1, p2 = case
    result = model_contact(1., 1., u1[0], u2[0], u1[1], u2[1],
                           p1[0], p2[0], p1[1], p2[1], m1, m2, k, mode)
    vx1, vx2, vy1, vy2, x1, x2, y1, y2 = result
    return (vx1, vy1), (vx2, vy2), (x1, y1), (x2, y2)


def number(value, spec=".6f"):
    mantissa, separator, exponent = format(value, spec).partition("e")
    if "." in mantissa:
        mantissa = mantissa.rstrip("0").rstrip(".")
    if mantissa == "-0":
        mantissa = "0"
    return mantissa + separator + exponent


def vector(v):
    return f"({number(v[0])}; {number(v[1])})"


def main():
    lines = ["# Результаты экспериментов M2", "",
             "Воспроизведение: `python docs/experiments.py`. Все величины в СИ. Радиусы обоих шаров — 1 м, расчёт начинается в момент касания. Коэффициент k = 10⁸: Н/м для Гука, Н/м³ᐟ² для Герца.", ""]
    max_energy_error, max_momentum_error = 0., 0.
    for mode, name in [(1, "Гук"), (2, "Герц")]:
        lines += [f"## Закон: {name}", "",
                  "| Опыт | m₁; m₂, кг | u₁; u₂, м/с | v₁ расчёт / теория, м/с | v₂ расчёт / теория, м/с | Относительная ошибка E |",
                  "|---|---|---|---|---|---:|"]
        for case in CASES:
            title, m1, m2, u1, u2, p1, p2 = case
            v1, v2, _, _ = simulate_case(case, mode)
            ref1, ref2 = elastic_velocities(m1, m2, u1, u2, p1, p2)
            error = abs(energy(m1, m2, v1, v2) / energy(m1, m2, u1, u2) - 1.)
            p_error = max(abs(a-b) for a, b in zip(momentum(m1,m2,v1,v2), momentum(m1,m2,u1,u2)))
            max_energy_error = max(max_energy_error, error)
            max_momentum_error = max(max_momentum_error, p_error)
            lines.append(f"| {title} | {m1:g}; {m2:g} | {vector(u1)}; {vector(u2)} | {vector(v1)} / {vector(ref1)} | {vector(v2)} / {vector(ref2)} | {number(error, '.2e')} |")
        lines.append("")
    lines += [f"Максимальная относительная ошибка энергии: {number(max_energy_error, '.3e')}; максимальная абсолютная ошибка компоненты импульса: {number(max_momentum_error, '.3e')} кг·м/с.", "",
              "## Жёсткость и приближение мгновенного удара", "",
              "Косой удар равных масс: u₁ = (2; 1), u₂ = (0; 0) м/с. Теоретически v₁ = (0; 1), v₂ = (2; 0) м/с. Ошибка — максимальное отклонение компоненты скорости от этого приближения.", "",
              "| Закон | k в единицах закона | Ошибка скорости, м/с |",
              "|---|---:|---:|"]
    case = CASES[4]
    ref1, ref2 = elastic_velocities(*case[1:])
    for mode, name in [(1, "Гук"), (2, "Герц")]:
        for k in [1e6, 1e8, 1e10]:
            v1, v2, _, _ = simulate_case(case, mode, k)
            error = max(abs(a-b) for a, b in zip(v1+v2, ref1+ref2))
            lines.append(f"| {name} | {k:g} | {number(error, '.6e')} |")
    lines += ["", "Эта серия меняет физическую жёсткость, а не только шаг интегрирования, поэтому не является отдельной проверкой порядка численного метода.", ""]
    target = Path(__file__).with_name("results.md")
    target.write_text("\n".join(lines), encoding="utf-8")
    print(f"Записан {target}")


if __name__ == "__main__":
    main()
