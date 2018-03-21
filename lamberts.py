import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler

from poliastro.iod import izzo

plt.rc('text', usetex=True)

x = np.linspace(-1, 2, num=1000)
M_list = 0, 1, 2, 3
ll_list = 1, 0.9, 0.7, 0, -0.7, -0.9, -1

fig, ax = plt.subplots(figsize=(8,6))
ax.set_prop_cycle(cycler('linestyle', ['-', '--']) * (cycler('color', ['black']) * len(ll_list)))

for M in M_list:
    for ll in ll_list:
        T_x0 = np.zeros_like(x)
        for ii in range(len(x)):
            y = izzo._compute_y(x[ii],ll)
            T_x0[ii] = izzo._tof_equation(x[ii], y, 0.0, ll, M)
        if M == 0 and ll == 1:
            T_x0[x > 0]=np.nan
            l, =ax.plot(x, T_x0)

            fig, ax = plt.subplots(figsize=(8, 6))
ax.set_prop_cycle(cycler('linestyle', ['-', '--']) *
                  (cycler('color', ['black']) * len(ll_list)))
for M in M_list:
    for ll in ll_list:
        T_x0 = np.zeros_like(x)
        for ii in range(len(x)):
            y = izzo._compute_y(x[ii], ll)
            T_x0[ii] = izzo._tof_equation(x[ii], y, 0.0, ll, M)
        if M == 0 and ll == 1:
            T_x0[x > 0] = np.nan
        elif M > 0:
            # Mask meaningless solutions
            T_x0[x > 1] = np.nan
        l, = ax.plot(x, T_x0)

ax.set_ylim(0, 10)

ax.set_xticks((-1, 0, 1, 2))
ax.set_yticks((0, np.pi, 2 * np.pi, 3 * np.pi))
ax.set_yticklabels(('$0$', '$\pi$', '$2 \pi$', '$3 \pi$'))

ax.vlines(1, 0, 10)
ax.text(0.65, 4.0, "elliptic")
ax.text(1.16, 4.0, "hyperbolic")

ax.text(0.05, 1.5, "$M = 0$", bbox=dict(facecolor='white'))
ax.text(0.05, 5, "$M = 1$", bbox=dict(facecolor='white'))
ax.text(0.05, 8, "$M = 2$", bbox=dict(facecolor='white'))

ax.annotate("$\lambda = 1$", xy=(-0.3, 1), xytext=(-0.75, 0.25), arrowprops=dict(arrowstyle="simple", facecolor="black"))
ax.annotate("$\lambda = -1$", xy=(0.3, 2.5), xytext=(0.65, 2.75), arrowprops=dict(arrowstyle="simple", facecolor="black"))

ax.grid()
ax.set_xlabel("$x$")
ax.set_ylabel("$T$")

for M in M_list:
    for ll in ll_list:
        x_T_min, T_min = izzo._compute_T_min(ll, M, 10, 1e-8)
        ax.plot(x_T_min, T_min, 'kx', mew=2)

T_ref = 1
ll_ref = 0

(x_ref, _), = izzo._find_xy(ll_ref, T_ref, 0, 10, 1e-8)
print(x_ref)

ax.plot(x_ref, T_ref, 'o', mew=2, mec='red', mfc='none')

from astropy import units as u
from poliastro.bodies import Earth
k = Earth.k
r0 = [15945.34, 0.0, 0.0] * u.km
r = [12214.83399, 10249.46731, 0.0] * u.km
tof = 76.0 * u.min

expected_va = [2.058925, 2.915956, 0.0] * u.km / u.s
expected_vb = [-3.451569, 0.910301, 0.0] * u.km / u.s

(v0, v), = izzo.lambert(k, r0, r, tof)
print(v)
