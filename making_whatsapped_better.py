import scipy.integrate as scint
import numpy as np
import matplotlib.pyplot as plt

# parameters for model: initial state, coefficients, timeframe
# (unreduced (p, q)-stationary unidirectional simple model)
u0_input = [0.5, 0.5, 1, 1, 1, 1]  # p, q, x, y, r, s
coeffs_input = [0.2, 1, 1, 0.2]  # k1, k2, k3, k4
time_bound = 50 # time until which to keep solving the ODE
time_resolution = 50 # number of points in time to actually log the values at (not counting t=0)
custom_timestamps = False # whether to use our own fixed logging points or let the solver decide for itself

def f(t, u, coeffs):  # differentiaalvergelijkingen
    p, q, x, y, r, s = u
    k1, k2, k3, k4 = coeffs
    A = np.array([
        0,
        0,
        k1 * p - k2 * q * x + k3 * x**2 * y - k4 * x,
        k2 * q * x - k3 * x**2 * y,
        k2 * q * x,
        k4 * x
    ])
    return(A)


# solving the ODE
time_span_input = (0, time_bound)
timestamps_input = np.linspace(0, time_bound, time_resolution + 1) if custom_timestamps else None
result = scint.solve_ivp(fun=f,
                         t_span=time_span_input,
                         y0=u0_input,
                         args=(coeffs_input,),
                         method='Radau',
                         t_eval=timestamps_input)
print(result)

# plotting
t = result.t
print(f"{len(t)=}")
p, q, x, y, r, s = result.y

plt.plot(t, p)
plt.plot(t, q)
plt.plot(t, x)
plt.plot(t, y)
plt.plot(t, r)
plt.plot(t, s)
plt.show()
