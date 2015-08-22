import sympy as sp
import numpy as np
import matplotlib.pylab as plt

sp.init_printing()

from scipy.integrate import odeint as sp_odeint


def f(y, t):
    return [y[1], -np.sin(y[0]) - np.sin(y[0] - 2 * t)]


t0 = 0
t1 = 200
dt = .1
y0 = [-.5, 0]

t2 = np.linspace(0, 200, 20001)
y2 = sp_odeint(f, y0, t2)
plt.plot(t2, y2)
plt.show()
