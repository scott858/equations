import sympy as sp
import numpy as np
import matplotlib.pylab as plt

x, v = sp.symbols('x v')
f = sp.Function('f')(x, v)
expr = -x - v
f = sp.lambdify((x, v), expr)

x = np.linspace(-1, 1, 10)
v = np.linspace(-1, 1, 10)
x, v = np.meshgrid(x, v)

result = f(x, v)
scale = np.sqrt(1 / (v ** 2 + result ** 2))

plt.quiver(x, v, v * scale, result * scale)
plt.show()

