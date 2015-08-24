import sympy as sp
import numpy as np
import matplotlib.pylab as plt

sp.init_printing()

from scipy.integrate import odeint as sp_odeint


def v_(z, t):
    return [z[1], -np.sin(z[0]) - np.sin(z[0] - 2 * t)]


N = 40
ti = 0
tf = 50
x0 = -.5
v0 = 0
T = np.linspace(ti, tf, (tf - ti) + 1)

delta_zi = 10 ** -5 * (np.random.rand(N, 2) - .5)
delta_zi[0, :] = [0, 0]

z0_bar = delta_zi.copy()
z0_bar[:, 0] += x0
z0_bar[:, 1] += v0

normed_delta_zi = np.square(delta_zi)
normed_delta_zi = np.sum(normed_delta_zi, 1)
normed_delta_zi = np.sqrt(normed_delta_zi)

z_lambda_list = list()

for time_interval in T[1:]:
    zf_list = list()
    t = np.linspace(0, time_interval, 2001)
    for z0 in z0_bar:
        zf_list.append(sp_odeint(v_, z0, t)[-1])

    zf = np.array(zf_list)
    normed_zf = zf.copy()
    normed_zf[:, 0] = normed_zf[:, 0] - normed_zf[0, 0]
    normed_zf[:, 1] = normed_zf[:, 1] - normed_zf[0, 1]
    normed_zf = np.square(normed_zf)
    normed_zf = np.sum(normed_zf, 1)
    normed_zf = np.sqrt(normed_zf)

    z_lambda_list.append(np.sum(np.log(normed_zf[1:] / normed_delta_zi[1:])) / ((N - 1) * time_interval))

z_lambda = np.array(z_lambda_list)
plt.plot(z_lambda)
plt.show()
