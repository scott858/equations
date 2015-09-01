import pyximport; pyximport.install()
from integrate import integrate, sin_cy2
from math import pi

def main():
    a, b = 0.0, 2.0 * pi
    return integrate(a, b, sin_cy2, N=400000)

if __name__ == '__main__':
    import cProfile
    cProfile.run('main()', sort='time')
