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

from euler import RungeKutta

sp.init_printing()

from scipy.integrate import odeint as sp_odeint

sys.setrecursionlimit(1000000)
resource.setrlimit(resource.RLIMIT_STACK, (2 ** 29, -1))

ti = 0
tf = 2000
time_interval = tf - ti
dt = .1
G = 10 ** -0

number_particles = 4
number_dimensions = 3
phase_index = 2
particle_index = number_particles * number_dimensions


def r(i, z):
    start_index = i * number_dimensions
    end_index = (i + 1) * number_dimensions
    r = z[0, 0, start_index:end_index]
    return r


def v(i, z):
    start_index = i * number_dimensions
    end_index = (i + 1) * number_dimensions
    v = z[0, 1, start_index:end_index]
    return v


def central_force(i, j, z):
    # exerted on i by j
    displacement = r(j, z) - r(i, z)
    mag_displacement = np.dot(displacement, displacement)
    mag_displacement = np.sqrt(mag_displacement)
    force = G * displacement / (mag_displacement ** (3 / 2))
    return force


def f(t, z):
    v = z[0, 1, :].copy()
    force = np.zeros_like(v)
    for i in range(number_particles):
        start_index = i * number_dimensions
        end_index = (i + 1) * number_dimensions
        particle_force = np.zeros(number_dimensions)
        for j in range(0, i):
            particle_force += central_force(i, j, z)
        for j in range(i + 1, number_particles):
            particle_force += central_force(i, j, z)
        force[start_index:end_index] = particle_force

    return np.array([[v, force]])


z0 = np.array([[[-10, 0, 0, 10, 0, 0, 0, 10, 0, 0, -10, 0], [0, -1, 0, 0, 1, 0, 1, 0, 0, -1, 0, 0]]])

rk = RungeKutta(delta_t=dt, time_interval=time_interval, w0=z0)
rk.set_f(f)
rk.z()
positions = rk.z_array[:, 0, :]

x1 = rk.z_array[:, 0, 0]
y1 = rk.z_array[:, 0, 1]
z1 = rk.z_array[:, 0, 2]

x2 = rk.z_array[:, 0, 3]
y2 = rk.z_array[:, 0, 4]
z2 = rk.z_array[:, 0, 5]

x3 = rk.z_array[:, 0, 6]
y3 = rk.z_array[:, 0, 7]
z3 = rk.z_array[:, 0, 7]

x4 = rk.z_array[:, 0, 8]
y4 = rk.z_array[:, 0, 9]
z4 = rk.z_array[:, 0, 10]

x = np.array([x1, x2, x3, x4])
y = np.array([y1, y2, y3, y4])
z = np.array([z1, z2, z3, z4])

trajectories = []
fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1], projection='3d')

colors = plt.cm.jet(np.linspace(0, 1, number_particles))
lines = sum([ax.plot([], [], [], '.', c=c) for c in colors], [])

lim = 20
ax.set_xlim3d([-lim, lim])
ax.set_ylim3d([-lim, lim])
ax.set_zlim3d([-lim, lim])
# time_template = 'time = %.1fs'
# time_text = ax.text(0.05, 0.9, 0.05, '', transform=ax.transAxes)


def init():
    for line in lines:
        line.set_data([], [])
        line.set_3d_properties([])

    return lines


stride = 5
def animate(i):
    xi = x[:, :i]
    yi = y[:, :i]
    zi = z[:, :i]
    for j, line in enumerate(lines, start=0):
        line.set_data(xi[j, -1:], yi[j, -1:])
        line.set_3d_properties(zi[j, -1:])

    fig.canvas.draw()
    # time_text.set_text(time_template % (i * dt))
    return lines


trajectory_index = np.arange(1, len(x1), 5)
ani = animation.FuncAnimation(fig, animate, trajectory_index, interval=1, blit=True, init_func=init)
# ani.save('particles.mp4', writer='mencoder', fps=15)
plt.show()
