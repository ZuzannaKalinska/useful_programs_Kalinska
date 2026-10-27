import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from pathlib import Path

# Parametry (Indeks 374061)
v0 = 37.0  # m/s
alpha_deg = 61.0  # stopnie
alpha = np.radians(alpha_deg)
g = 9.81  # m/s^2
m = 1.0   # kg
k = 0.1   # współczynnik oporu powietrza (kg/s dla liniowego lub kg/m dla kwadratowego)

# Katalog na pliki
out_dir = Path("pages/solutions/04_sol_vscode/files/rzut")
out_dir.mkdir(parents=True, exist_ok=True)

# 1. Rzut bez oporu powietrza (analitycznie)
t_max_ideal = 2 * v0 * np.sin(alpha) / g
t_ideal = np.linspace(0, t_max_ideal, 100)
x_ideal = v0 * np.cos(alpha) * t_ideal
y_ideal = v0 * np.sin(alpha) * t_ideal - 0.5 * g * t_ideal**2

# 2. Rzut z oporem powietrza (numerycznie - model kwadratowy F = -k*v*|v|)
def projectile_with_drag(t, state, k, m, g):
    x, y, vx, vy = state
    v = np.sqrt(vx**2 + vy**2)
    ax = -(k/m) * v * vx
    ay = -g - (k/m) * v * vy
    return [vx, vy, ax, ay]

v0x = v0 * np.cos(alpha)
v0y = v0 * np.sin(alpha)
state0 = [0, 0, v0x, v0y]

# Rozwiązanie do momentu uderzenia w ziemię (y=0)
def hit_ground(t, state, *args):
    return state[1]
hit_ground.terminal = True
hit_ground.direction = -1

sol = solve_ivp(projectile_with_drag, [0, 10], state0, args=(k, m, g), 
                events=hit_ground, t_eval=np.linspace(0, 10, 500))

x_drag = sol.y[0]
y_drag = sol.y[1]

# Wykres 1: Porównanie trajektorii
plt.figure(figsize=(10, 6))
plt.plot(x_ideal, y_ideal, 'r--', label='Bez oporu powietrza')
plt.plot(x_drag, y_drag, 'b-', label='Z oporem powietrza (kwadratowym)')
plt.xlabel('Dystans x [m]')
plt.ylabel('Wysokość y [m]')
plt.title(f'Trajektoria rzutu ukośnego ($v_0={v0}$ m/s, $\\alpha={alpha_deg}^\\circ$)')
plt.legend()
plt.grid(True)
plt.savefig(out_dir / "trajektoria.png", dpi=150)
plt.close()

# Wykres 2: Prędkość w czasie
plt.figure(figsize=(10, 6))
v_ideal = np.sqrt((v0 * np.cos(alpha))**2 + (v0 * np.sin(alpha) - g * t_ideal)**2)
vx_drag = sol.y[2]
vy_drag = sol.y[3]
v_drag = np.sqrt(vx_drag**2 + vy_drag**2)

plt.plot(t_ideal, v_ideal, 'r--', label='v (bez oporu)')
plt.plot(sol.t, v_drag, 'b-', label='v (z oporem)')
plt.xlabel('Czas t [s]')
plt.ylabel('Prędkość v [m/s]')
plt.title('Zmiana prędkości w czasie')
plt.legend()
plt.grid(True)
plt.savefig(out_dir / "predkosc.png", dpi=150)
plt.close()

# Wykres 3: Energia
plt.figure(figsize=(10, 6))
Ek_ideal = 0.5 * m * v_ideal**2
Ep_ideal = m * g * y_ideal
Etot_ideal = Ek_ideal + Ep_ideal

Ek_drag = 0.5 * m * v_drag**2
Ep_drag = m * g * y_drag
Etot_drag = Ek_drag + Ep_drag

plt.plot(sol.t, Etot_drag, 'b-', label='Energia całkowita (z oporem)')
plt.plot(t_ideal, Etot_ideal, 'r--', label='Energia całkowita (bez oporu)')
plt.xlabel('Czas t [s]')
plt.ylabel('Energia [J]')
plt.title('Zachowanie energii całkowitej')
plt.legend()
plt.grid(True)
plt.savefig(out_dir / "energia.png", dpi=150)
plt.close()

print("Wykresy zostały wygenerowane.")
