from cython cimport boundscheck, wraparound

def summer(double[:] mv):
    cdef:
        int i, N
        double ss = 0.0

    N = mv.shape[0]
    with boundscheck(False), wraparound(False):
        for i in range(N):
            ss += mv[i]
    return ss