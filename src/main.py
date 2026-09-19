
from math import *
import matplotlib.pyplot as plt

def time_to_contact(r1, r2, v1x, v2x, v1y, v2y, x1, x2, y1, y2):
    s = r1 + r2;
    ux = (v1x - v2x)
    uy = (v1y - v2y)

    dx = (x1 - x2)
    dy = (y1 - y2)

    A = ux ** 2 + uy ** 2
    B = (ux * dx + uy * dy)
    C = dx ** 2 + dy **2 - s ** 2;

    D = B ** 2 - A * C;

    if (A == 0 or B>= 0 or D < 0 or C < 0):
        return None;
    return C / (-B + D**0.5)


def model_contact(r1, r2, v1x, v2x, v1y, v2y, x1, x2, y1, y2, m1, m2, k, mode):
    Fx = 0
    Fy = 0
    dt = sqrt(((m1 * m2) / (m1 + m2) )/ k) / 10000.0
    T = [0.0]
    n = 0
    while True:
        delta = (r1 + r2) - dist((x1, y1), (x2, y2))
        deltax_norm = (x1 - x2) / hypot(x1- x2, y1 - y2)
        deltay_norm = (y1 - y2) / hypot(x1- x2, y1 - y2)
        if (delta < 0):
            delta = 0
        if (mode == 1):
            Fx = - delta * k * deltax_norm
            Fy = - delta * k * deltay_norm
        if (mode == 2):
            Fx = - (delta ** (3/2)) * k * deltax_norm
            Fy = - (delta ** (3/2)) * k * deltay_norm
        v1x -= (Fx * dt) / m1
        v1y -= (Fy * dt) / m1
        v2x += (Fx * dt) / m2
        v2y += (Fy * dt) / m2
        x1 += v1x * dt
        y1 += v1y * dt
        x2 += v2x * dt
        y2 += v2y * dt
        if (delta <= 0 and n > 0):
            break
        n += 1;
    return v1x, v2x, v1y, v2y, x1, x2, y1, y2


def draw(r1, r2, v1x, v2x, v1y, v2y, start1, start2, touch1, touch2, x_contact, y_contact):
    t_after = 1.5
    end1 = (touch1[0] + v1x * t_after, touch1[1] + v1y * t_after)
    end2 = (touch2[0] + v2x * t_after, touch2[1] + v2y * t_after)

    figure, axis = plt.subplots(figsize=(10, 6))

    axis.plot([start1[0], touch1[0]], [start1[1], touch1[1]], "C0")
    axis.plot([start2[0], touch2[0]], [start2[1], touch2[1]], "C1")
    axis.plot([touch1[0], end1[0]], [touch1[1], end1[1]], "C0--")
    axis.plot([touch2[0], end2[0]], [touch2[1], end2[1]], "C1--")
    axis.plot(x_contact, y_contact, "ko")
    for c, r, color, style in [(start1, r1, "C0", "-"), (start2, r2, "C1", "-"),
                               (touch1, r1, "C0", "--"), (touch2, r2, "C1", "--"),
                               (end1, r1, "C0", ":"), (end2, r2, "C1", ":")]:
        axis.add_patch(plt.Circle(c, r, fill=False, color=color, linestyle=style, alpha=0.6))

    axis.set_xlabel("x")
    axis.set_ylabel("y")
    axis.set_aspect("equal")
    axis.grid(alpha=0.3)
    axis.legend()
    plt.show()

def main():
    # Шар 1
    x1, y1 = -30.0, 0.0
    v1x, v1y = 20.0, 0.0
    r1 = 5.0
    m1 = 4.0
    # Шар 2
    x2, y2 = 19.0, 10.0
    v2x, v2y = -10.0, -10.0
    r2 = 5.0
    m2 = 5.0

    start1, start2 = (x1, y1), (x2, y2)

    k = 1e8          # жёсткость
    mode = 1         # 1 - закон Гука, 2 - закон Герца

    t = time_to_contact(r1, r2, v1x, v2x, v1y, v2y, x1, x2, y1, y2)
    if t is None:
        print("Столкновение не произошло")
        return
    x1 += t * v1x
    y1 += t * v1y
    x2 += t * v2x
    y2 += t * v2y
    touch1, touch2 = (x1, y1), (x2, y2)
    x_contact = x1 + (x2 - x1) * (r1 / (r1 + r2))
    y_contact = y1 + (y2 - y1) * (r1 / (r1 + r2))

    v1x, v2x, v1y, v2y, x1, x2, y1, y2 = model_contact(r1, r2, v1x, v2x, v1y, v2y, x1, x2, y1, y2, m1, m2, k, mode)
    print(f"Столкновение произошло через {t:.2f} секунд")
    print(f"Координата столкновения: ({x_contact:.2f}, {y_contact:.2f})")
    print(f"Первый шар: V1 = ({v1x:.2f}, {v1y:.2f})")
    print(f"Второй шар: V2 = ({v2x:.2f}, {v2y:.2f})")
    draw(r1, r2, v1x, v2x, v1y, v2y, start1, start2, touch1, touch2, x_contact, y_contact)



if __name__ == "__main__":
    main()
