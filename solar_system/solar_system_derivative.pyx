# cython: profile=True
from libc.math cimport pow, sqrt
from solar_system.solar_system_data import *

phase_index = 2
start_v_index = number_particles * number_dimensions


cdef double[:] r_cy(int i, double[:] w):
    cdef:
        int start_index, end_index
    start_index = i * number_dimensions
    end_index = (i + 1) * number_dimensions

    cdef double[:] ri = np.empty(end_index - start_index, dtype=np.double)
    ri = w[start_index:end_index]
    return ri


cdef double[:] my_subtract(double[:] x, double[:] y):
    cdef:
        int i
        double[:] x_minus_y = np.empty(x.shape[0], dtype=np.double)

    for i in range(x.shape[0]):
        x_minus_y[i] = x[i] - y[i]

    return x_minus_y


cdef double[:] central_force(int i, int j, double[:] w):
    cdef:
        double mj, mag_displacement = 0
        double force_coefficient
        int k
    # exerted on i by j
    rj = r_cy(j, w)
    ri = r_cy(i, w)
    displacement = my_subtract(rj, ri)

    for k in range(displacement.shape[0]):
        mag_displacement += displacement[k] * displacement[k]
    mag_displacement = sqrt(mag_displacement)

    mj = masses[j]
    force_coefficient = G * mj / (pow(mag_displacement, 3))

    for k in range(displacement.shape[0]):
        displacement[k] = force_coefficient * displacement[k]
    return displacement


cdef double[:] f_core(double[:] w):
    cdef:
        int v_start, i, j, start_index, end_index
        double[:] particle_force = np.empty(number_dimensions, dtype=np.double)
        double[:] central_force_ij = np.empty(number_dimensions, dtype=np.double)
        double[:] w_ = np.empty(w.shape[0], np.double)
    v_start = number_dimensions * number_particles
    w_[:v_start] = w[v_start:]
    for i in range(number_particles):
        particle_force[0] = 0
        particle_force[1] = 0
        particle_force[2] = 0

        start_index = i * number_dimensions + v_start
        end_index = (i + 1) * number_dimensions + v_start
        for j in range(0, i):
            central_force_ij = central_force(i, j, w)
            particle_force[0] += central_force_ij[0]
            particle_force[1] += central_force_ij[1]
            particle_force[2] += central_force_ij[2]
        for j in range(i + 1, number_particles):
            central_force_ij = central_force(i, j, w)
            particle_force[0] += central_force_ij[0]
            particle_force[1] += central_force_ij[1]
            particle_force[2] += central_force_ij[2]
        w_[start_index:end_index] = particle_force
    return w_


def f(w, t):
    # f = dv / dt
    return np.asarray(f_core(np.asarray(w, dtype=np.double)), dtype=np.double)


def v(int i, w):
    cdef:
        int start_index, end_index
    start_index = i * number_dimensions
    end_index = (i + 1) * number_dimensions
    vi = w[start_v_index + start_index:end_index]
    return vi


def r(int i, w):
    cdef:
        int start_index, end_index
    start_index = i * number_dimensions
    end_index = (i + 1) * number_dimensions

    ri = w[start_index:end_index]
    return ri


