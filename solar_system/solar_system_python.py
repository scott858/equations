import sympy as sp
import numpy as np
import matplotlib

matplotlib.use('QT4Agg')
import matplotlib.pyplot as plt
import sys
import resource
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from .solar_system_data import *

sp.init_printing()

from scipy.integrate import odeint as sp_odeint

# sys.setrecursionlimit(1000000)
# resource.setrlimit(resource.RLIMIT_STACK, (2 ** 29, -1))

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

number_particles = len(masses)
number_dimensions = 3
phase_index = 2
start_v_index = number_particles * number_dimensions


def r(i, w):
    start_index = i * number_dimensions
    end_index = (i + 1) * number_dimensions
    ri = w[start_index:end_index]
    return ri


def v(i, w):
    start_index = i * number_dimensions
    end_index = (i + 1) * number_dimensions
    vi = w[start_v_index + start_index:end_index]
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
    # f = dv / dt
    v_start = number_dimensions * number_particles
    v_slice = w[v_start:]
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


w0 = np.array([], dtype=np.float64)
for x0i, y0i, z0i in zip(x0, y0, z0):
    w0 = np.concatenate([w0, [x0i], [y0i], [z0i]])
for vx0i, vy0i, vz0i in zip(vx0, vy0, vz0):
    w0 = np.concatenate([w0, [vx0i], [vy0i], [vz0i]])

t = np.arange(0, time_interval, dt, dtype=np.float64)


def run():
    # w, s = sp_odeint(f, w0, t, full_output=True)
    w = sp_odeint(f, w0, t)
    # print(w[-1])
    # print(s)

    # fig = plt.figure()
    # ax = fig.add_axes([0, 0, 1, 1], projection='3d')
    # ax.patch.set_facecolor('black')
    # ax.axis('off')
    #
    # colors = (
    # '#ffff00', '#cc0000', '#ff0000', '#00ffff', '#cc0066', '#663300', '#ff3300', '#ff6699', '#9933ff', '#33ccff',)
    # widths = np.log(1 + masses / m_sun)
    # widths = widths / np.min(widths)
    # widths = 1 + np.log(widths)
    # widths = np.floor(widths)
    # lines = sum([ax.plot([], [], [], '.', c=c, ms=width) for c, width in zip(colors, widths)], [])
    #
    # lim = 20 * earth_distance
    #
    # ax.set_xlim3d([-lim, lim])
    # ax.set_ylim3d([-lim, lim])
    # # lim = lim / 10
    # ax.set_zlim3d([-lim, lim])
    # time_template = 'time = %.1fdays'
    # time_text = ax.text(0.05, 0.9, 0.05, '', transform=ax.transAxes)
    # time_text.set_color('#ffffff')
    #
    # def init():
    #     for line in lines:
    #         line.set_data([], [])
    #         line.set_3d_properties([])
    #     return lines
    #
    # def animate(i):
    #     for j, line in enumerate(lines, start=0):
    #         x = w[i, number_dimensions * j]
    #         y = w[i, number_dimensions * j + 1]
    #         z = w[i, number_dimensions * j + 2]
    #         line.set_data(x, y)
    #         line.set_3d_properties(z)
    #
    #     fig.canvas.draw()
    #     time_text.set_text(time_template % (i * dt / (seconds_per_hour * hours_per_day)))
    #     return lines

    # stride = 5
    # trajectory_index = np.arange(1, w.shape[0], 10)
    # ani = animation.FuncAnimation(fig, animate, trajectory_index, interval=1, blit=True)
    # ani.save('particles.mp4', writer='mencoder', fps=15)

    # for c, j in zip(colors, range(number_particles)):
    #     # all
    #     x = w[:, number_dimensions * j]
    #     y = w[:, number_dimensions * j + 1]
    #     z = w[:, number_dimensions * j + 2]
    #     ax.plot(x, y, z, c)

    # plt.show()
    return w
