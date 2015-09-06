# cython: profile=True
include 'solar_system_data_cy.pxi'
# from solar_system.solar_system_data import *
from libc.math cimport pow, sqrt


cpdef get_t():
    return np.arange(0, time_interval, dt, dtype=np.double)


cpdef get_w0():
    w0 = np.array([], dtype=np.double)
    for x0i, y0i, z0i in zip(x0, y0, z0):
        w0 = np.concatenate([w0, [x0i], [y0i], [z0i]])
    for vx0i, vy0i, vz0i in zip(vx0, vy0, vz0):
        w0 = np.concatenate([w0, [vx0i], [vy0i], [vz0i]])
    return w0


cdef double[:] r(int i, double[:] w):
    cdef:
        int start_index, end_index
    start_index = i * number_dimensions
    end_index = (i + 1) * number_dimensions
    return w[start_index:end_index]


cdef double[:] displacement = np.empty(number_dimensions, dtype=np.double)
cdef double[:] central_force(int i, int j, double[:] w):
    cdef:
        double mj, mag_displacement = 0
        double force_coefficient
        int k
    # exerted on i by j
    rj = r(j, w)
    ri = r(i, w)
    for k in range(rj.shape[0]):
        displacement[k] = rj[k] - ri[k]

    for k in range(displacement.shape[0]):
        mag_displacement += displacement[k] * displacement[k]
    mag_displacement = sqrt(mag_displacement)

    mj = masses[j]
    force_coefficient = G * mj / (pow(mag_displacement, 3))

    for k in range(displacement.shape[0]):
        displacement[k] = force_coefficient * displacement[k]
    return displacement

cdef:
    double[:] particle_force = np.empty(number_dimensions, dtype=np.double)
    double[:] central_force_ij = np.empty(number_dimensions, dtype=np.double)
    double[:] w_ = np.empty(2 * number_dimensions * number_particles, np.double)
    int v_start = number_dimensions * number_particles
cdef double[:] f_core(double[:] w):
    cdef:
        int i, j, start_index, end_index
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
