from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp


# Parametry wynikajace z numeru indeksu 374061.
v0 = 37.0
alpha = np.deg2rad(61.0)
g = 9.81
m = 1.0
k = 0.1  # kg/m; F_d = -k * |v| * v

base_dir = Path(__file__).parent
out_dir = base_dir / "files" / "rzut"
out_dir.mkdir(parents=True, exist_ok=True)

v0x = v0 * np.cos(alpha)
v0y = v0 * np.sin(alpha)

# Ruch idealny: wyniki analityczne.
t_flight_ideal = 2 * v0y / g
t_ideal = np.linspace(0, t_flight_ideal, 600)
x_ideal = v0x * t_ideal
y_ideal = v0y * t_ideal - 0.5 * g * t_ideal**2
vx_ideal = np.full_like(t_ideal, v0x)
vy_ideal = v0y - g * t_ideal
speed_ideal = np.hypot(vx_ideal, vy_ideal)


def projectile_with_drag(time, state):
    _, _, vx, vy = state
    _ = time
    speed = np.hypot(vx, vy)
    ax = -(k / m) * speed * vx
    ay = -g - (k / m) * speed * vy
    return vx, vy, ax, ay


def hit_ground(time, state):
    _ = time
    return state[1]


hit_ground.terminal = True
hit_ground.direction = -1

solution = solve_ivp(
    projectile_with_drag,
    (0, 20),
    (0, 0, v0x, v0y),
    events=hit_ground,
    dense_output=True,
    rtol=1e-9,
    atol=1e-11,
    max_step=0.02,
)

if not solution.success or len(solution.t_events[0]) == 0:
    raise RuntimeError("Nie udalo sie wyznaczyc chwili upadku ciala.")

t_flight_drag = solution.t_events[0][0]
t_drag = np.linspace(0, t_flight_drag, 600)
x_drag, y_drag, vx_drag, vy_drag = solution.sol(t_drag)
speed_drag = np.hypot(vx_drag, vy_drag)
ax_drag = -(k / m) * speed_drag * vx_drag
ay_drag = -g - (k / m) * speed_drag * vy_drag

def save_plot(filename, title, xlabel, ylabel):
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / filename, dpi=180)
    plt.savefig(out_dir / f"{Path(filename).stem}.svg", format="svg")
    plt.close()


# 1. Tor ruchu i punkty charakterystyczne.
plt.figure(figsize=(9, 5.5))
plt.plot(x_ideal, y_ideal, "--", label="Bez oporu")
plt.plot(x_drag, y_drag, label="Z oporem kwadratowym")
plt.scatter([x_ideal.max(), x_drag.max()], [0, 0], zorder=3, label="Punkt upadku")
save_plot("trajektoria.png", "Tor rzutu ukośnego", "Położenie x [m]", "Wysokość y [m]")

# 2. Współrzędne w czasie.
plt.figure(figsize=(9, 5.5))
plt.plot(t_ideal, x_ideal, "--", label="x, bez oporu")
plt.plot(t_drag, x_drag, label="x, z oporem")
plt.plot(t_ideal, y_ideal, ":", label="y, bez oporu")
plt.plot(t_drag, y_drag, label="y, z oporem")
save_plot("polozenie_czas.png", "Położenie w funkcji czasu", "Czas t [s]", "Położenie [m]")

# 3. Składowe prędkości.
plt.figure(figsize=(9, 5.5))
plt.plot(t_ideal, vx_ideal, "--", label="$v_x$, bez oporu")
plt.plot(t_drag, vx_drag, label="$v_x$, z oporem")
plt.plot(t_ideal, vy_ideal, ":", label="$v_y$, bez oporu")
plt.plot(t_drag, vy_drag, label="$v_y$, z oporem")
save_plot("skladowe_predkosci.png", "Składowe prędkości", "Czas t [s]", "Prędkość [m/s]")

# 4. Moduł prędkości.
plt.figure(figsize=(9, 5.5))
plt.plot(t_ideal, speed_ideal, "--", label="Bez oporu")
plt.plot(t_drag, speed_drag, label="Z oporem")
save_plot("szybkosc.png", "Moduł prędkości w czasie", "Czas t [s]", "Szybkość |v| [m/s]")

# 5. Przyspieszenie w modelu z oporem.
plt.figure(figsize=(9, 5.5))
plt.plot(t_drag, ax_drag, label="$a_x$, z oporem")
plt.plot(t_drag, ay_drag, label="$a_y$, z oporem")
plt.axhline(-g, color="black", linestyle="--", label="$-g$, bez oporu")
save_plot("przyspieszenie.png", "Składowe przyspieszenia", "Czas t [s]", "Przyspieszenie [m/s²]")

# 6. Energia kinetyczna, potencjalna i mechaniczna.
kinetic_ideal = 0.5 * m * speed_ideal**2
potential_ideal = m * g * y_ideal
mechanical_ideal = kinetic_ideal + potential_ideal
kinetic_drag = 0.5 * m * speed_drag**2
potential_drag = m * g * y_drag
mechanical_drag = kinetic_drag + potential_drag

plt.figure(figsize=(9, 5.5))
plt.plot(t_ideal, mechanical_ideal, "--", label="Mechaniczna, bez oporu")
plt.plot(t_drag, mechanical_drag, label="Mechaniczna, z oporem")
plt.plot(t_drag, kinetic_drag, ":", label="Kinetyczna, z oporem")
plt.plot(t_drag, potential_drag, "-.", label="Potencjalna, z oporem")
save_plot("energia.png", "Energie w czasie", "Czas t [s]", "Energia [J]")

ideal_range = x_ideal[-1]
drag_range = x_drag[-1]
ideal_height = y_ideal.max()
drag_height = y_drag.max()
print(f"Model idealny: lot {t_flight_ideal:.3f} s, zasieg {ideal_range:.3f} m, wysokosc {ideal_height:.3f} m")
print(f"Model z oporem: lot {t_flight_drag:.3f} s, zasieg {drag_range:.3f} m, wysokosc {drag_height:.3f} m")
print(f"Wykresy zapisano w: {out_dir}")
