import sympy as sp
import numpy as np
import matplotlib

matplotlib.use('QT4Agg')
import matplotlib.pyplot as plt
import sys
import resource
from matplotlib.colors import cnames
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from solar_system_data import *

from euler import RungeKutta

sp.init_printing()

from scipy.integrate import odeint as sp_odeint

sys.setrecursionlimit(1000000)
resource.setrlimit(resource.RLIMIT_STACK, (2 ** 29, -1))

ti = 0
years = 1
days_per_year = 365
hours_per_day = 24
seconds_per_hour = 3600
tf = years * days_per_year * hours_per_day * seconds_per_hour * time_scale
time_interval = tf - ti
dt = 24 * 3600 * time_scale
G = 6.67 * 10 ** -11
G = (space_scale * space_scale) * G / (mass_scale * time_scale ** 2)

number_particles = 10
number_dimensions = 3
phase_index = 2
particle_index = number_particles * number_dimensions


def r(i, w):
    start_index = i * number_dimensions
    end_index = (i + 1) * number_dimensions
    r = w[start_index:end_index].copy()
    return r


def v(i, w):
    start_index = i * number_dimensions
    end_index = (i + 1) * number_dimensions
    vi = w[particle_index + start_index:end_index].copy()
    return vi


def central_force(i, j, w):
    # exerted on i by j
    displacement = r(j, w) - r(i, w)
    mag_displacement = np.dot(displacement, displacement)
    mag_displacement = np.sqrt(mag_displacement)
    mj = masses[j]
    v_ = G * displacement * mj / (mag_displacement ** 3)
    return v_


def f(w, t):
    v_start = number_dimensions * number_particles
    v_slice = w[v_start:].copy()
    v_ = np.zeros_like(v_slice, dtype=np.float64)
    for i in range(number_particles):
        start_index = i * number_dimensions
        end_index = (i + 1) * number_dimensions
        particle_force = np.zeros(number_dimensions, dtype=np.float64)
        for j in range(0, i):
            particle_force += central_force(i, j, w)
        for j in range(i + 1, number_particles):
            particle_force += central_force(i, j, w)
        v_[start_index:end_index] = particle_force
    return np.concatenate([v_slice, v_])

# wx0 = 50 * (np.random.rand(particle_index) - .5)
# wv0 = 50 * (np.random.rand(particle_index) - .5)
# w0 = np.concatenate([wx0, wv0])
w0 = np.array([], dtype=np.float64)
for x0i, y0i, z0i in zip(x0, y0, z0):
    w0 = np.concatenate([w0, [x0i], [y0i], [z0i]])
for vx0i, vy0i, vz0i in zip(vx0, vy0, vz0):
    w0 = np.concatenate([w0, [vx0i], [vy0i], [vz0i]])

# w0 = np.array([x0[0], y0[0], z0[0], x0[3], y0[3], z0[3], vx0[0], vy0[0], vz0[0], vx0[3], vy0[3], vz0[3]], dtype=np.float64)

# w0 = np.concatenate([[-10, 0, 0, 10, 0, 0, 0, 10, 0, 0, -10, 0], [0, -1, 0, 0, 1, 0, 1, 0, 0, -1, 0, 0]])
# w0 = np.concatenate([[-10, 0, 0, 10, 0, 0], [0, -1, 0, 0, 1, 0]])
t = np.arange(0, time_interval, dt, dtype=np.float64)
w_result = sp_odeint(f, w0, t, mxstep=50000000, hmin=.0001)

trajectories = []
fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1], projection='3d')
# ax.axis('off')

colors = plt.cm.jet(np.linspace(0, 1, number_particles))
lines = sum([ax.plot([], [], [], '.', c=c) for c in colors], [])

lim = 2.0 * earth_distance
ax.set_xlim3d([-lim, lim])
ax.set_ylim3d([-lim, lim])
ax.set_zlim3d([-lim, lim])
time_template = 'time = %.1fdays'
time_text = ax.text(0.05, 0.9, 0.05, '', transform=ax.transAxes)


def init():
    for line in lines:
        line.set_data([], [])
        line.set_3d_properties([])

    return lines


stride = 5


def animate(i):
    for j, line in enumerate(lines, start=0):
        x = w_result[i, number_dimensions * j]
        y = w_result[i, number_dimensions * j + 1]
        z = w_result[i, number_dimensions * j + 2]
        line.set_data(x, y)
        line.set_3d_properties(z)

    fig.canvas.draw()
    time_text.set_text(time_template % (i * dt / (seconds_per_hour * hours_per_day)))
    return lines


trajectory_index = np.arange(1, w_result.shape[0], 1)
ani = animation.FuncAnimation(fig, animate, trajectory_index, interval=1, blit=True)
# ani.save('particles.mp4', writer='mencoder', fps=15)
plt.show()
