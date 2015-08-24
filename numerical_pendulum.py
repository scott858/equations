import sympy as sp
import numpy as np
import matplotlib.pylab as plt

sp.init_printing()

from scipy.integrate import odeint as sp_odeint


def v_(z, t):
    return np.array([z[1], -np.sin(z[0]) - np.sin(z[0] - 2 * t)])


ti = 0
tf = 200
dt = .1
xi_vi = [-.5, 0]

t2 = np.linspace(ti, tf, 20001)
t1 = np.linspace(ti, tf, 2001)

x1_v1 = sp_odeint(v_, xi_vi, t1)
x2_v2 = sp_odeint(v_, xi_vi, t2)

plt.plot(t2, x2_v2)
plt.plot(t1, x1_v1)
plt.show()
