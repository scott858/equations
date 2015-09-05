# cython: profile=True
import sympy as sp
import numpy as np
import matplotlib

matplotlib.use('QT4Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

from solar_system.solar_system_data import *
from solar_system.solar_system_derivative import f

sp.init_printing()

from scipy.integrate import odeint as sp_odeint

w0 = np.array([], dtype=np.double)
for x0i, y0i, z0i in zip(x0, y0, z0):
    w0 = np.concatenate([w0, [x0i], [y0i], [z0i]])
for vx0i, vy0i, vz0i in zip(vx0, vy0, vz0):
    w0 = np.concatenate([w0, [vx0i], [vy0i], [vz0i]])

t = np.arange(0, time_interval, dt, dtype=np.double)


def run():
    # w, s = sp_odeint(f, w0, t, full_output=True)
    w = sp_odeint(f, w0, t)
    return w


def solar_animate(w):
    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1], projection='3d')
    ax.patch.set_facecolor('black')
    ax.axis('off')

    colors = (
        '#ffff00', '#cc0000', '#ff0000', '#00ffff', '#cc0066', '#663300', '#ff3300', '#ff6699', '#9933ff', '#33ccff',)
    widths = np.log(1 + masses / m_sun)
    widths = widths / np.min(widths)
    widths = 1 + np.log(widths)
    widths = np.floor(widths)
    lines = sum([ax.plot([], [], [], '.', c=c, ms=width) for c, width in zip(colors, widths)], [])

    lim = 20 * earth_distance

    ax.set_xlim3d([-lim, lim])
    ax.set_ylim3d([-lim, lim])
    # lim = lim / 10
    ax.set_zlim3d([-lim, lim])
    time_template = 'time = %.1fdays'
    time_text = ax.text(0.05, 0.9, 0.05, '', transform=ax.transAxes)
    time_text.set_color('#ffffff')

    def init():
        for line in lines:
            line.set_data([], [])
            line.set_3d_properties([])
        return lines

    def animate(i):
        for j, line in enumerate(lines, start=0):
            x = w[i, number_dimensions * j]
            y = w[i, number_dimensions * j + 1]
            z = w[i, number_dimensions * j + 2]
            line.set_data(x, y)
            line.set_3d_properties(z)

        fig.canvas.draw()
        time_text.set_text(time_template % (i * dt / (seconds_per_hour * hours_per_day)))
        return lines

    trajectory_index = np.arange(1, w.shape[0], 10)
    ani = animation.FuncAnimation(fig, animate, trajectory_index, interval=1, blit=True)
    ani.save('particles.mp4', writer='mencoder', fps=15)

    plt.show()


def solar_plot(w):
    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1], projection='3d')
    ax.patch.set_facecolor('black')
    ax.axis('off')

    colors = (
        '#ffff00', '#cc0000', '#ff0000', '#00ffff', '#cc0066', '#663300', '#ff3300', '#ff6699', '#9933ff', '#33ccff',)

    for c, j in zip(colors, range(number_particles)):
        # all
        x = w[:, number_dimensions * j]
        y = w[:, number_dimensions * j + 1]
        z = w[:, number_dimensions * j + 2]
        ax.plot(x, y, z, c)


if __name__ == '__main__':
    import cProfile
    cProfile.run('run()', sort='time')

