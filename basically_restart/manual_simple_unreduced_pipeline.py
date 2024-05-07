from new_simulator import simulate
from visualize_simple_unreduced_1 import visualize_simple_unreduced_1
import numpy as np


def f(t, u, coeffs):
    p, q, x, y, r, s = u
    k1, k2, k3, k4 = coeffs
    return np.array([
        0,
        0,
        k1 * p - k2 * q * x + k3 * x * x * y - k4 * x,
        k2 * q * x - k3 * x * x * y,
        k2 * q * x,
        k4 * x
    ])


initials = [0.59, 0.56, 0.83, 1.25, 0, 0]  # list of starting values of the variables; first part of the parameters.
coefficients = [0.19, 1.07, 0.85, 0.22]  # list of coefficients for reaction speeds; second part of the parameters.
interval = (0, 200)  # cutoff point in time to stop the simulation at, or None for the default value of 50.
granularity = None  # number of points in time to actually log the values at (not counting t=0),
# or None to let the solver itself decide for us.

if granularity is None:
    evaluations_list = None
else:
    evaluations_list = [np.linspace(interval[0], interval[1], granularity + 1)]

result = simulate(f, [initials], [coefficients], [interval], evaluations_list)
f, initials_list, coefficients_list, interval_list, evaluations_list, solution_list = result

visualize_simple_unreduced_1(f, initials_list, coefficients_list, interval_list, evaluations_list, solution_list)
