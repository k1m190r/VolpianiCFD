import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc

rc("font", family="serif", size=16)
rc("lines", linewidth=1.5)
plt.rc("legend", **{"fontsize": 14})

Nx = 101
xmax = 2.0
xmin = -2.0
Lx = xmax - xmin
dx = Lx / (Nx - 1)
x = np.linspace(xmin, xmax, Nx)

dt = 0.04
t_end = 5.0

Nt = int(t_end / dt)
a = 0.8
CFL = a * dt / dx
U = np.zeros((Nt + 1, Nx))
U[0, :] = np.exp(-0.5 * (x / 0.4) ** 2)
Uex = U[0, :]

for n in range(0, Nt):
    if a > 0.0:
        for i in range(1, Nx):
            U[n + 1, i] = U[n, i] - CFL * (U[n, i] - U[n, i - 1])
        U[n + 1, 0] = U[n + 1, Nx - 1]
    else:
        for i in range(0, Nx - 1):
            U[n + 1, i] = U[n, 1] - CFL * (U[n, i + 1] - U[n, i])
        U[n + 1, Nx - 1] = U[n, 0]

    d = a * (n + 1) * dt
    Uex = np.exp(-0.5 * (np.mod(x - d + xmax, 4) - xmax) ** 2 / 0.4**2)

    if n == 0:
        fig, ax = plt.subplots(figsize=(5.5, 4))
    plt.clf()
    plt.plot(x, U[n + 1, :])
