cdef class Particle:

    cdef readonly double mass, position, velocity

    def __init__(self, m, p, v):
        self.mass = m
        self.position = p
        self.velocity = v

    def get_momentum(self):
        return self.mass * self.velocity

    @staticmethod
    def add_momenta(particles):
        total_mom = 0.0
        for particle in particles:
            total_mom += particle.get_momentum()

        return total_mom

    cdef double get_momentum_c(self):
        return self.mass * self.velocity

    @staticmethod
    def add_momenta_c(list particles):
        cdef:
            double total_mom = 0.0
            Particle particle

        for particle in particles:
            total_mom += particle.get_momentum_c()

        return total_mom
