import sympy as sp
import numpy as np
import matplotlib.pylab as plt

sp.init_printing()

x = sp.Function('x')
t = sp.Symbol('t')
w0 = sp.Symbol('w0', real=True, positive=True)
C0, C1, C2, C3 = sp.symbols('C0 C1 C2 C3')

x1 = sp.Function('x1')
x2 = sp.Function('x2')
eq1 = sp.Eq(x1(t).diff(t, t), -x1(t) + 2 * (x2(t) - x1(t)))
eq2 = sp.Eq(x2(t).diff(t, t), -x2(t) + 2 * (x1(t) - x2(t)))
s = sp.dsolve([eq1, eq2])
x1 = s[0].rhs
x2 = s[1].rhs
print(x1)
print(x2)
ic1 = sp.Eq(x1.subs(t, 0), 0)
ic2 = sp.Eq(x1.diff(t).subs(t, 0), 0)
ic3 = sp.Eq(x2.subs(t, 0), 1)
ic4 = sp.Eq(x2.diff(t).subs(t, 0), 0)

C = sp.Symbol('C')
C = sp.solve([ic1, ic2, ic3, ic4])
print(C)

x1c = x1.subs(C)
x1c = x1c.simplify()
x1cl = sp.lambdify(t, x1c, 'numpy')

x2c = x2.subs(C)
x2c = x2c.simplify()
x2cl = sp.lambdify(t, x2c, "numpy")

npt = np.linspace(0, 100, 1000)
plt.plot(x1cl(npt))
plt.plot(x2cl(npt))
plt.show()
