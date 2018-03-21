import matplotlib.pyplot as plt
plt.ion()
from astropy import units as u
from astropy.time import Time
from astropy.utils.data import conf
conf.remote_timeout = 10000
from astropy.coordinates import solar_system_ephemeris
solar_system_ephemeris.set("jpl")
from poliastro.bodies import *
from poliastro.twobody import Orbit
from poliastro.plotting import OrbitPlotter, plot

EPOCH = Time("2017-09-01 12:05:50", scale="tdb")
earth = Orbit.from_body_ephem(Earth, EPOCH)
plot(earth, label=Earth)

