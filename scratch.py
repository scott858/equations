import sympy as sp
import numpy as np
import matplotlib.pylab as plt
import sys
import resource

from euler import RungeKutta

sp.init_printing()

from scipy.integrate import odeint as sp_odeint

sys.setrecursionlimit(1000000)
resource.setrlimit(resource.RLIMIT_STACK, (2 ** 29, -1))

ti = 0
tf = 200
time_interval = tf - ti
dt = .2

number_particles = 4
number_dimensions = 3
time_index = 1
phase_index = 2
particle_index = number_particles * number_dimensions

z0 = np.random.rand(time_index, phase_index, particle_index)


def r(i, n, z):
    start_index = i * number_dimensions
    end_index = (i + 1) * number_dimensions
    return z[n, 0, start_index:end_index]


def v(i, n, z):
    start_index = i * number_dimensions
    end_index = (i + 1) * number_dimensions
    return z[n, 1, start_index:end_index]


def central_force(i, j, n, z):
    # exerted on i by j
    displacement = r(j, n, z) - r(i, n, z)
    force = displacement / ((displacement.T * displacement) ** (3 / 2))
    return force


def f(n, z):
    v = z[n, 1, :].copy()
    force = np.zeros_like(v)
    for i in range(number_particles):
        start_index = i * number_dimensions
        end_index = (i + 1) * number_dimensions
        particle_force = np.zeros(number_dimensions)
        for j in range(0, i):
            particle_force += central_force(i, j, n, z)
        for j in range(i + 1, number_particles):
            particle_force += central_force(i, j, n, z)
        force[start_index:end_index] = particle_force

    return [z[n, 1, :], force]

print(z0)
print('\n')
print(r(0, 0, z0))
print('\n')
print(v(0, 0, z0))
print('\n')
print(central_force(0, 1, 0, z0))
print('\n')
print(f(0, z0))

# rk = RungeKutta(delta_t=dt, time_interval=time_interval, z0=z0)
# rk.set_f(chaotic_pendulum)
#
# rk.v()
# plt.plot(rk.time_array, rk.v_array)
# plt.show()
