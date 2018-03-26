# %matplotlib inline
from matplotlib.axes import Axes
from matplotlib.pyplot import gcf, gca
import matplotlib.pyplot as plt
from matplotlib.lines import *
from matplotlib import animation, rc
import numpy as np
from astropy import units as u
from poliastro.plotting import plot, OrbitPlotter
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.maneuver import Maneuver

rc('animation', html='html5')
    
op = OrbitPlotter()
ss_i = Orbit.circular(Earth, alt=700 * u.km)
hoh = Maneuver.hohmann(ss_i, 36000 * u.km)
ss_a, ss_f = ss_i.apply_maneuver(hoh, intermediate=True)
op.plot(ss_i, label="Initial orbit")
op.plot(ss_a, label="Transfer orbit")
ln = op.plot(ss_f, label="Final orbit")[0]
ln
xlim = op.ax.get_xlim()
ylim = op.ax.get_xlim()
xdata, ydata = [], []

fig = gcf()

def project(op, rr):
    rr_proj = rr - rr.dot(op._frame[2])[:, None] * op._frame[2]
    x = rr_proj.dot(op._frame[0])
    y = rr_proj.dot(op._frame[1])
    return x, y

rr = ss_f.sample().get_xyz().transpose()
x, y = project(op, rr)

# initialization function: plot the background of each frame
def init():
    op.ax.set_xlim(xlim)
    op.ax.set_ylim(ylim)
    return ln

# animation function. This is called sequentially
def animate(i):
    # a, = op.ax.plot(x.to(u.km).value[i], y.to(u.km).value[i], '.')
    xdata.append(x.to(u.km).value[i])
    ydata.append(y.to(u.km).value[i])
    ln.set_data(xdata, ydata)
    return ln

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=10, interval=20, blit=True)

plt.show()
