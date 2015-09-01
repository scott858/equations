cimport cython

from libc.math cimport sin as sin_cy
from math import sin

def integrate_py(double a, double b, f, int N=2000):
    dx = (b - a) / N
    s = 0
    for i in range(N):
        s += f(a + i * dx)
    return s * dx


def sin_py2(x):
    return sin(x)**2


def sin_cy2(x):
    return sin_cy(x)**2


@cython.cdivision(True)
cdef double integrate_cy(double a, double b, double (*f)(double), int N=2000):

    cdef:
        int i
        double dx = (b - a) / N
        double s = 0.0

    for i in range(N):
        s += f(a + i * dx)
    return s * dx
