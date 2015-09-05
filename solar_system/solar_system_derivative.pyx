from .solar_system_data import *

phase_index = 2
start_v_index = number_particles * number_dimensions


cdef r(int i, w):
    cdef:
        int start_index, end_index
    start_index = i * number_dimensions
    end_index = (i + 1) * number_dimensions
    ri = w[start_index:end_index]
    return ri


cdef v(int i, w):
    cdef:
        int start_index, end_index
    start_index = i * number_dimensions
    end_index = (i + 1) * number_dimensions
    vi = w[start_v_index + start_index:end_index]
    return vi


cdef central_force(int i, int j, w):
    cdef:
        double mag_displacement, mj
    # exerted on i by j
    displacement = r(j, w) - r(i, w)
    mag_displacement = np.dot(displacement, displacement)
    mag_displacement = np.sqrt(mag_displacement)
    mj = masses[j]
    v_ = G * displacement * mj / (mag_displacement ** 3)
    return v_


cpdef f(w, t):
    # f = dv / dt
    cdef:
        int v_start, i, j
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


