import numpy as np
import matplotlib.pyplot as plt
plt.ion()

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris

from poliastro.bodies import Sun, Earth, Jupiter
from poliastro.twobody import Orbit
from poliastro.maneuver import Maneuver
from poliastro.iod import izzo
from poliastro.plotting import plot, OrbitPlotter
from poliastro.util import norm
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

solar_system_ephemeris.set("jpl")
## Initial data
# Links and sources: https://github.com/poliastro/poliastro/wiki/EuroPython:-Per-Python-ad-Astra
date_launch = Time("2011-08-05 16:25", scale='utc')
C_3 = 31.1 * u.km**2 / u.s**2
date_flyby = Time("2013-10-09 19:21", scale='utc')
date_arrival = Time("2016-07-05 03:18", scale='utc')
# Initial state of the Earth
ss_e0 = Orbit.from_body_ephem(Earth, date_launch)
r_e0, v_e0 = ss_e0.rv()

# State of the Earth the day of the flyby
ss_efly = Orbit.from_body_ephem(Earth, date_flyby)
r_efly, v_efly = ss_efly.rv()

# Assume that the insertion velocity is tangential to that of the Earth
dv = C_3**.5 * v_e0 / norm(v_e0)
man = Maneuver.impulse(dv)
# Inner Cruise 1
ic1 = ss_e0.apply_maneuver(man)
ic1.rv()
# We propagate until the aphelion
ss_aph = ic1.propagate(ic1.period / 2)
# Let's compute the Lambert solution to do the flyby of the Earth
time_of_flight = date_flyby - ss_aph.epoch
(v_aph, v_fly), = izzo.lambert(Sun.k, ss_aph.r, ss_efly.r, time_of_flight)
# Check the delta-V
norm(v_aph - ss_aph.v)  # Too high!

ss_aph_post = Orbit.from_vectors(Sun, ss_aph.r, v_aph, epoch=ss_aph.epoch)
ss_junofly = Orbit.from_vectors(Sun, r_efly, v_fly, epoch=date_flyby)


ss_j = Orbit.from_body_ephem(Jupiter, date_arrival)
r_j, v_j = ss_j.rv()

(v_flypre, v_oip), = izzo.lambert(Sun.k, r_efly, r_j, date_arrival - date_flyby)
ss_oip = Orbit.from_vectors(Sun, r_j, v_oip, epoch=date_flyby)
fig, ax = plt.subplots(figsize=(9, 12))

op = OrbitPlotter(ax)

op.plot(ss_e0, label="Earth")
op.plot(ic1, label="Inner Cruise 1")
#op.plot(ss_efly)
op.plot(ss_aph_post, label="Back to Earth")
op.plot(ss_oip, label="Jupiter Orbit Insertion Phase")
op.plot(ss_j, label="Jupiter")

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
rr_e0 = ss_e0.sample().get_xyz().transpose().to(u.km).value
rr_ic1,_ = split_list(ic1.sample().get_xyz().transpose().to(u.km).value) # Only get half because of Hohmann Maneuver
rr_aph_post = ss_aph_post.sample().get_xyz().transpose().to(u.km).value
rr_oip,_ = split_list(ss_oip.sample().get_xyz().transpose().to(u.km).value) # Only get half
rr_j = ss_j.sample().get_xyz().transpose().to(u.km).value
rr = np.concatenate((rr_e0, rr_ic1, rr_aph_post, rr_oip, rr_j))
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
ani = FuncAnimation(fig, update, frames=np.linspace(0, len(x) - 1, len(x)), init_func=init, blit=True)
plt.show()
