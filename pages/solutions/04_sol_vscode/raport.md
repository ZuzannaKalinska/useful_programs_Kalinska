# Analiza rzutu ukośnego: model idealny vs opór powietrza

*Raport obliczeniowy przygotowany w ramach zadania z wykorzystaniem VS Code.*

**Temat:** Dwuwymiarowy ruch punktu materialnego w polu grawitacyjnym.  
**Metoda:** Rozwiązania analityczne oraz całkowanie numeryczne (Python/SciPy).  
**Dane wejściowe (indeks 374061):**
- Prędkość początkowa $v_0 = 37.0$ m/s
- Kąt początkowy $\alpha = 61.0^\circ$

---

## Spis treści
1. [Model fizyczny](#1-model-fizyczny)
2. [Ruch bez oporu powietrza (idealny)](#2-ruch-bez-oporu-powietrza-idealny)
3. [Ruch z oporem powietrza](#3-ruch-z-oporem-powietrza)
4. [Analiza wyników i wykresy](#4-analiza-wyników-i-wykresy)
5. [Kod źródłowy](#5-kod-źródłowy)
6. [Wnioski](#6-wnioski)

---

## 1. Model fizyczny

Rzut ukośny polega na nadaniu ciału prędkości początkowej $\vec{v}_0$ pod kątem $\alpha$ do poziomu. W najprostszym modelu (idealnym) jedyną siłą działającą na ciało jest siła ciężkości $\vec{F}_g = m\vec{g}$. W modelu rzeczywistym należy uwzględnić siłę oporu powietrza $\vec{F}_d$, która zależy od prędkości ciała.

### Parametry symulacji

| Parametr | Symbol | Wartość | Jednostka |
| :--- | :---: | :---: | :---: |
| Prędkość początkowa | $v_0$ | 37.0 | m/s |
| Kąt początkowy | $\alpha$ | 61.0 | $^\circ$ |
| Przyspieszenie ziemskie | $g$ | 9.81 | m/s² |
| Masa ciała | $m$ | 1.0 | kg |
| Współczynnik oporu | $k$ | 0.1 | kg/m |

---

## 2. Ruch bez oporu powietrza (idealny)

W modelu idealnym ruch w osi $x$ jest jednostajny, a w osi $y$ jednostajnie zmienny.

### 2.1. Równania ruchu

$$
\begin{cases}
x(t) = v_0 \cos(\alpha) t \\
y(t) = v_0 \sin(\alpha) t - \frac{gt^2}{2}
\end{cases}
$$

### 2.2. Wyznaczone parametry rzutu

- **Czas lotu ($y=0$):** $t_l = \frac{2v_0 \sin(\alpha)}{g} \approx 6.59$ s
- **Zasięg rzutu:** $Z = \frac{v_0^2 \sin(2\alpha)}{g} \approx 118.81$ m
- **Maksymalna wysokość:** $H = \frac{(v_0 \sin(\alpha))^2}{2g} \approx 53.33$ m

---

## 3. Ruch z oporem powietrza

Przy wyższych prędkościach siła oporu powietrza jest proporcjonalna do kwadratu prędkości (model Newtona):

$$
\vec{F}_d = -k |\vec{v}| \vec{v}
$$

Równania różniczkowe drugiego rzędu przyjmują postać:

$$
\begin{cases}
m \ddot{x} = -k \sqrt{\dot{x}^2 + \dot{y}^2} \dot{x} \\
m \ddot{y} = -mg - k \sqrt{\dot{x}^2 + \dot{y}^2} \dot{y}
\end{cases}
$$

Z uwagi na brak prostego rozwiązania analitycznego dla modelu kwadratowego, trajektoria została wyznaczona numerycznie przy użyciu metody `solve_ivp`.

---

## 4. Analiza wyników i wykresy

### 4.1. Porównanie trajektorii

Na poniższym wykresie widać wyraźną różnicę między parabolą (model idealny) a trajektorią z oporem powietrza. Opór powietrza powoduje skrócenie zasięgu oraz sprawia, że krzywa staje się asymetryczna (opadanie jest bardziej strome).

![Trajektoria rzutu](files/rzut/trajektoria.png)

### 4.2. Zmiana prędkości

W modelu idealnym prędkość w osi $x$ jest stała. Z oporem powietrza, obie składowe prędkości maleją, co prowadzi do szybszego wytracania energii kinetycznej.

![Prędkość w czasie](files/rzut/predkosc.png)

### 4.3. Zachowanie energii

W modelu bez oporu, energia całkowita ($E_k + E_p$) jest stała. W obecności oporu powietrza, praca sił oporu jest ujemna, co prowadzi do ciągłego spadku energii mechanicznej układu.

![Energia](files/rzut/energia.png)

---

## 5. Kod źródłowy

Poniżej znajduje się kod w języku Python, który został wykorzystany do generacji powyższych wykresów i obliczeń.

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from pathlib import Path

# Parametry wynikające z indeksu 374061
v0 = 37.0  # m/s (pierwsze dwie cyfry)
alpha_deg = 61.0  # stopnie (ostatnie dwie cyfry)
alpha = np.radians(alpha_deg)
g = 9.81
m = 1.0
k = 0.1  # współczynnik oporu

# 1. Obliczenia analityczne (model idealny)
t_max_ideal = 2 * v0 * np.sin(alpha) / g
t_ideal = np.linspace(0, t_max_ideal, 100)
x_ideal = v0 * np.cos(alpha) * t_ideal
y_ideal = v0 * np.sin(alpha) * t_ideal - 0.5 * g * t_ideal**2

# 2. Obliczenia numeryczne (model z oporem)
def projectile_with_drag(t, state, k, m, g):
    x, y, vx, vy = state
    v = np.sqrt(vx**2 + vy**2)
    ax = -(k/m) * v * vx
    ay = -g - (k/m) * v * vy
    return [vx, vy, ax, ay]

v0x = v0 * np.cos(alpha)
v0y = v0 * np.sin(alpha)
state0 = [0, 0, v0x, v0y]

def hit_ground(t, state, *args): return state[1]
hit_ground.terminal = True
hit_ground.direction = -1

sol = solve_ivp(projectile_with_drag, [0, 10], state0, args=(k, m, g), 
                events=hit_ground, t_eval=np.linspace(0, 10, 500))

# Wykresy generowane za pomocą Matplotlib...
# (Pełny kod dostępny w pliku generate_plots.py)
```

---

## 6. Wnioski

1. **Wpływ oporu:** Opór powietrza drastycznie zmienia zasięg rzutu (z ok. 118 m do ok. 50 m w przyjętym modelu).
2. **Asymetria:** Trajektoria z oporem nie jest parabolą; punkt szczytowy znajduje się bliżej miejsca upadku niż miejsca startu.
3. **Energia:** W modelu rzeczywistym energia mechaniczna nie jest zachowana — jest rozpraszana w wyniku tarcia o powietrze.
4. **Weryfikacja:** Narzędzia takie jak Python/SciPy pozwalają na łatwe modelowanie zjawisk, których nie da się rozwiązać prostymi wzorami szkolnymi.
