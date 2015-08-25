import sympy as sp
import numpy as np
import matplotlib.pylab as plt
import sys
import resource
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

from euler import RungeKutta

sp.init_printing()

from scipy.integrate import odeint as sp_odeint

sys.setrecursionlimit(1000000)
resource.setrlimit(resource.RLIMIT_STACK, (2 ** 29, -1))

ti = 0
tf = 200
time_interval = tf - ti
dt = .01
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

z0 = np.random.rand(1, phase_index, particle_index) - .5
z0 = np.array([[[-10, 0, 0, 10, 0, 0, 0, 10, 0, 0, -10, 0], [0, -1, 0, 0, 1, 0, 1, 0, 0, -1, 0, 0]]])
# z0 = np.array([[[-100, 0, 100, 0, 0, 100, 0, -100], [0, -1, 0, 1, 1, 0, -1, 0]]])

rk = RungeKutta(delta_t=dt, time_interval=time_interval, z0=z0)
rk.set_f(f)
rk.z()
plt.subplots(1)
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

fig2 = plt.figure()
ax = fig2.add_subplot(111, autoscale_on=False, projection='3d')
# ax = Axes3D(fig2)
ax.grid()
lim = 100
ax.set_xlim3d([-lim, lim])
ax.set_ylim3d([-lim, lim])
ax.set_zlim3d([-lim, lim])

line, = ax.plot([], [], [], '.')
time_template = 'time = %.1fs'
time_text = ax.text(0.5, 0.9, 0.5, '', transform=ax.transAxes)


def init():
    line.set_data(x1[0], y1[0])
    line.set_3d_properties(z1[0])
    time_text.set_text('')
    return line, time_text


def animate(i):
    x = [x1[:i], x2[:i], x3[:i], x4[:i]]
    y = [y1[:i], y2[:i], y3[:i], y4[:i]]
    z = [z1[:i], z2[:i], z3[:i], z4[:i]]

    line.set_data(z1, y1)
    line.set_3d_properties(z1)
    time_text.set_text(time_template % (i * dt))
    return line, time_text

# ani = animation.FuncAnimation(fig2, animate, np.arange(1, len(y1)), interval=.1, blit=True)
ani = animation.FuncAnimation(fig2, animate, np.arange(1, len(x1)), interval=.1, blit=True, init_func=init)

# ani.save('particles.mp4', writer='mencoder', fps=15)
plt.show()

# print(z0)
# print('\n')
# print(r(0, 0, z0))
# print('\n')
# print(v(0, 0, z0))
# print('\n')
# print(central_force(0, 1, 0, z0))
# print('\n')
# print(f(0, z0))

# rk = RungeKutta(delta_t=dt, time_interval=time_interval, z0=z0)
# rk.set_f(f)
# rk.z()

# print(rk.time_array)
# print('\n')
# print(rk.z_array)
# print('\n')
# print(rk.z_array[:, 0, 1])
# print('\n')
# print(rk.z_array.shape)
# plt.subplots(1)
# plt.plot(rk.time_array, rk.z_array[:, 0, :])
# plt.subplots(1)
# plt.plot(rk.time_array, rk.z_array[:, 1, :])
# plt.show()



# def v_(t, z):
#     return np.array([z[0, 1], -np.sin(z[0, 0]) - np.sin(z[0, 0] - 2 * t)])
#
#
# ti = 0
# tf = 200
# time_interval = tf - ti
# dt = .1
# z0 = np.array([[[-.5], [0]]])
#
# rk = RungeKutta(delta_t=dt, time_interval=time_interval, z0=z0)
# rk.set_f(v_)
# rk.z()
# plt.plot(rk.time_array, rk.z_array[:, 0, :])
# plt.show()
