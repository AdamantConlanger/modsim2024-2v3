import scipy.integrate as sc
import numpy as np
import matplotlib.pyplot as plt

# beginvoorwaardes
p = 0.5
q = 0.5
k1 = 0.2
k2 = 1
k3 = 1
k4 = 0.2


def fun(t, param):  # differentiaalvergelijkingen
    x = param[2]
    y = param[3]
    f1 = 0
    f2 = 0
    f3 = k1 * p - k2 * q * x + k3 * x**2 * y - k4 * x
    f4 = k2 * q * x - k3 * x**2 * y
    f5 = k2 * q * x
    f6 = k4 * x
    A = np.array([f1, f2, f3, f4, f5, f6])
    return (A)


# oplossen vergelijkingen
x_pts = np.linspace(0, 50, 51)
result = sc.solve_ivp(fun, (0, 50), [0.5, 0.5, 1, 1, 1, 1], t_eval=x_pts)
print(result)

# plotten
t = result.t
p = result.y[0]
q = result.y[1]
x = result.y[2]
y = result.y[3]
r = result.y[4]
s = result.y[5]

plt.plot(t, p)
plt.plot(t, q)
plt.plot(t, x)
plt.plot(t, y)
plt.plot(t, r)
plt.plot(t, s)
plt.show()
