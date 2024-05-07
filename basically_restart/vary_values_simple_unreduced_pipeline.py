from new_simulator import simulate
from visualize_simple_unreduced_1 import visualize_simple_unreduced_1
from cartesian_product import cartesian_product_for_variations
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


base_initials = [0.59, 0.56, 0.83, 1.25, 0, 0]  # list of starting values of the variables; first part of the parameters.
base_coefficients = [0.19, 1.07, 0.85, 0.22]  # list of coefficients for reaction speeds; second part of the parameters.
interval = (0, 200)  # cutoff point in time to stop the simulation at, or None for the default value of 50.
granularity = None  # number of points in time to actually log the values at (not counting t=0),
# or None to let the solver itself decide for us.

vary_simultaneously = False # whether to entrywise combine the variations (True) or Cartesian them (False)
multiplicative = False # whether to apply variations multiplicatively (True) or additively (False)
variations_initials = [None, np.linspace(0, 0.2, 3), None, np.linspace(0, 0.1, 2), [0], None]
variations_coefficients = [None, [0, 0.1], None, None]

tmp = np.array([int(multiplicative)])
variations_initials = [tmp if item is None else np.array(item) for item in variations_initials]
variations_coefficients = [tmp if item is None else np.array(item) for item in variations_coefficients]


if vary_simultaneously:
    variations_combined = variations_initials + variations_coefficients
    parallel_length = max([len(item) for item in variations_combined])
    for index in range(len(variations_combined)):
        if len(variations_combined[index]) == 1:
            variations_combined[index] = np.repeat(variations_combined[index], parallel_length)
        elif len(variations_combined[index]) != parallel_length:
            raise ValueError("Arrays aren't all of the same length or singletons")
    transposition = np.transpose(np.array(variations_combined))
    len_initials = len(variations_initials)
    len_both = len_initials + len(variations_coefficients)
    variations = [(item[0:len_initials], item[len_initials:len_both]) for item in transposition]
else:
    variations = cartesian_product_for_variations(variations_initials, variations_coefficients)

if multiplicative:
    initials_list = [np.array(base_initials) * item[0] for item in variations]
    coefficients_list = [np.array(base_coefficients) * item[1] for item in variations]
else:
    initials_list = [np.array(base_initials) + item[0] for item in variations]
    coefficients_list = [np.array(base_coefficients) + item[1] for item in variations]

interval_list = [interval for item in variations]

if granularity is None:
    evaluations_list = None
else:
    evaluations_list = [np.linspace(interval[0], interval[1], granularity + 1) for item in variations]

result = simulate(f, initials_list, coefficients_list, interval_list, evaluations_list)
f, initials_list, coefficients_list, interval_list, evaluations_list, solution_list = result

visualize_simple_unreduced_1(f, initials_list, coefficients_list, interval_list, evaluations_list, solution_list, suppress_legend=True, suppress_ghosts=True)
