from cpython.array cimport array
from cython cimport boundscheck, wraparound

def summer(double[:] mv):
    """Sums its argument's contents."""
    cdef double d, ss = 0.0
    for d in mv:
        ss += d
    return ss


@boundscheck(False)
@wraparound(False)
def cy_summer(double[:] mv):
    """Sums its argument's contents."""
    cdef:
        double d, ss = 0.0
        int i, N

    N = mv.shape[0]
    for i in range(N):
        ss += mv[i]
    return ss
