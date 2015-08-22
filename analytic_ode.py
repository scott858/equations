import sympy as sp

sp.init_printing()

x = sp.Function('x')
t = sp.Symbol('t')
w0 = sp.Symbol('w0', real=True, positive=True)
deq = sp.Eq(x(t).diff(t, t), -w0 ** 2 * x(t))
solution = sp.dsolve(deq)
C0, C1, C2, C3 = sp.symbols('C0 C1 C2 C3')
print(sp.solve(solution, C1))
print(solution.subs(t, 0))
x = solution.rhs
print(x)
print(sp.solve(sp.Eq(x.diff(t).subs(t, 0), 1), C1))
