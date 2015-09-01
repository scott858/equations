import pyximport

pyximport.install()
from integrate import integrate_py
from integrate import integrate_cy
from integrate import sin_cy2
from integrate import sin_py2
from math import pi


def main_cy():
    a, b = 0.0, 2.0 * pi
    return integrate_cy(a, b, sin_cy2, N=4000000)

def main_py():
    a, b = 0.0, 2.0 * pi
    return integrate_py(a, b, sin_py2, N=4000000)


if __name__ == '__main__':
    import cProfile

    cProfile.run('main_cy()', sort='time')
