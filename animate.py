# %matplotlib inline
from matplotlib.axes import Axes
from matplotlib.lines import *
from matplotlib.pyplot import gcf, gca
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from astropy import units as u
from poliastro.plotting import plot, OrbitPlotter
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.maneuver import Maneuver

# Plot result
op = OrbitPlotter()
ss_i = Orbit.circular(Earth, alt=700 * u.km)
hoh = Maneuver.hohmann(ss_i, 36000 * u.km)
ss_a, ss_f = ss_i.apply_maneuver(hoh, intermediate=True)
op.plot(ss_i, label="Initial orbit")
op.plot(ss_a, label="Initial orbit")
op.plot(ss_f, label="Final orbit")

# Get limits
xlim = op.ax.get_xlim()
ylim = op.ax.get_xlim()

# Set data
xdata, ydata = [], []

# Create figure
fig, ax = plt.subplots()
ln, = plt.plot([], [], 'r-', animated=True)
plt.xlabel('x (km)')
plt.ylabel('y (km)')

# Set a project function (this is the same function from _project in plotting.py)
def project(op, rr):
    rr_proj = rr - rr.dot(op._frame[2])[:, None] * op._frame[2]
    x = rr_proj.dot(op._frame[0])
    y = rr_proj.dot(op._frame[1])
    return x, y

# Based on https://stackoverflow.com/questions/752308/split-list-into-smaller-lists
def split_list(a_list):
    half = int(np.floor(len(a_list)/2))
    return a_list[:half], a_list[half:]

# Treat the data
rr_i = ss_i.sample().get_xyz().transpose().to(u.km).value
rr_a,_ = split_list(ss_a.sample().get_xyz().transpose().to(u.km).value) # Only get half because of Hohmann Maneuver
rr_f = ss_f.sample().get_xyz().transpose().to(u.km).value
rr = np.concatenate((rr_i, rr_a, rr_f))
x, y = project(op, rr)

# Init data
def init():
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    return ln,

# Set update plot function
def update(frame):
    xdata.append(x[int(frame)])
    ydata.append(y[int(frame)])
    ln.set_data(xdata, ydata)
    return ln,

# Create
ani = FuncAnimation(fig, update, frames=np.linspace(0, len(x) - 1, len(x)), init_func=init, blit=True, save_count=len(x))

# Save
ani.save('trajectory.mp4')

# Show
# plt.show()
