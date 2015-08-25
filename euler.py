import sympy as sp
import numpy as np
import logging

logger = logging.getLogger(__name__)

sp.init_printing()


class Euler:
    @staticmethod
    def _f_default(t, v):
        return t - v

    def __init__(self, delta_t, time_interval, z0):
        self._delta_t = delta_t
        self._time_interval = time_interval
        self._z0 = z0
        self._z_array = self._z0.copy()

        self._time_steps = int(self._time_interval / self._delta_t)

        self._f = Euler._f_default

        self._time_array = np.linspace(0, time_interval, time_interval / delta_t + 1)

    @property
    def time_steps(self):
        return int(self._time_steps)

    @property
    def z_array(self):
        return self._z_array

    @property
    def time_array(self):
        return self._time_array

    @time_steps.setter
    def time_steps(self, steps):
        self._time_steps = int(steps)

    def set_f(self, f):
        self._f = f

    def _t(self, n):
        return n * self._delta_t

    def z(self):
        return self._z_recursion(self.time_steps)

    def _z_recursion(self, n):
        if not type(n) == int:
            raise ValueError('n must be an integer')
        elif n < 0:
            raise ValueError('n must be greater than 0')
        elif n == 0:
            return self._z_array
        else:
            try:
                return self._z_array[n, :, :]
            except IndexError as err:
                logger.info(err)

            zn = self._step_calc(n)
            self._z_array = np.append(self._z_array, zn, axis=0)
            return zn

    def _step_calc(self, n):
        dt = self._delta_t
        tn_1 = self._t(n - 1)
        zn_1 = self._z_recursion(n - 1)

        zn = zn_1 + dt * self._f(tn_1, zn_1)
        return zn


class PredictorCorrector(Euler):
    def __init__(self, delta_t, time_interval, z0):
        super().__init__(delta_t=delta_t, time_interval=time_interval, z0=z0)

    def _step_calc(self, n):
        dt = self._delta_t
        tn_1 = self._t(n - 1)
        tn = self._t(n)

        zn_1 = self._z_recursion(n - 1)

        zn_euler = zn_1 + dt * self._f(tn_1, zn_1)
        zn_pc = zn_1 + (dt / 2) * (self._f(tn_1, zn_1) + self._f(tn, zn_euler))
        return zn_pc


class RungeKutta(Euler):
    def __init__(self, delta_t, time_interval, z0):
        super().__init__(delta_t=delta_t, time_interval=time_interval, z0=z0)

    def _step_calc(self, n):
        dt = self._delta_t
        tn_1 = self._t(n - 1)
        zn_1 = self._z_recursion(n - 1)

        k1 = self._f(tn_1, zn_1)
        k2 = self._f(tn_1 + dt / 2, zn_1 + dt * k1 / 2)
        k3 = self._f(tn_1 + dt / 2, zn_1 + dt * k2 / 2)
        k4 = self._f(tn_1 + dt, zn_1 + dt * k3)
        zn = zn_1 + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        return zn
