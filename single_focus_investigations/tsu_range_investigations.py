from simulate import simulate
from tsu_range_visualize import visualize
from cartesian_product import cartesian_product
import numpy as np


def f(t, u, coeffs):
    x, y = u
    k1, k2, k3, k4, p, q = coeffs
    return np.array([
        k1 * p - k2 * q * x + k3 * x * x * y - k4 * x,
        k2 * q * x - k3 * x * x * y
    ])


# list of starting values of the variables; first part of the parameters.
base_initials = [0.83, 1.25]
# list of coefficients for reaction speeds; second part of the parameters.
base_coefficients = [0.19, 1.07, 0.85, 0.22, 0.59, 0.56]
interval = (0, 1000)  # cutoff point in time to stop the simulation at, or None for the default value of 50.
granularity = 500  # number of points in time to actually log the values at (not counting t=0),
# or None to let the solver itself decide for us.
plotted_interval = None  # time span to actually plot, as closed interval. or None for full plot.

vary_simultaneously = False  # whether to entrywise combine the variations (True) or Cartesian them (False)
multiplicative = True  # whether to apply variations multiplicatively (True) or additively (False)
variations_initials = [None, None]
variations_coefficients = [np.linspace(1.3, 1.35, 20), None, None, None, None, None]
is_focus_initials = [False, False]
is_focus_coefficients = [True, False, False, False, False, False]

initials_length = len(base_initials)
coefficients_length = len(base_coefficients)
base_variables = base_initials + base_coefficients
variations_variables = variations_initials + variations_coefficients
is_focus = is_focus_initials + is_focus_coefficients

tmp = np.array([int(multiplicative)])
variations_variables = [tmp if item is None else np.array(item) for item in variations_variables]
minima_variations = [min(item) for item in variations_variables]
maxima_variations = [max(item) for item in variations_variables]

if vary_simultaneously:
    parallel_length = max([len(item) for item in variations_variables])
    for index in range(len(variations_variables)):
        if len(variations_variables[index]) == 1:
            variations_variables[index] = np.repeat(variations_variables[index], parallel_length)
        elif len(variations_variables[index]) != parallel_length:
            raise ValueError("Arrays aren't all of the same length or singletons")
    variations = np.transpose(np.array(variations_variables))
else:
    variations = cartesian_product(*variations_variables)

if multiplicative:
    variables_list = [np.array(base_variables) * item for item in variations]
    minima_variables = np.array(base_variables) * np.array(minima_variations)
    maxima_variables = np.array(base_variables) * np.array(maxima_variations)
    extrema_variables = [(minima_variables[index], maxima_variables[index]) for index in range(len(base_variables))]
else:
    variables_list = [np.array(base_variables) + item for item in variations]
    minima_variables = np.array(base_variables) + np.array(minima_variations)
    maxima_variables = np.array(base_variables) + np.array(maxima_variations)
    extrema_variables = [(minima_variables[index], maxima_variables[index]) for index in range(len(base_variables))]

initials_list = [variables[:initials_length] for variables in variables_list]
coefficients_list = [variables[initials_length:] for variables in variables_list]

interval_list = [interval for item in variations]

if granularity is None:
    evaluations_list = None
else:
    evaluations_list = [np.linspace(interval[0], interval[1], granularity + 1) for item in variations]

result = simulate(f, initials_list, coefficients_list, interval_list, evaluations_list)
f, initials_list, coefficients_list, interval_list, evaluations_list, solution_list = result

visualize(f, initials_list, coefficients_list, interval_list,
          evaluations_list, solution_list, is_focus, extrema_variables, suppress_legend=False, suppress_ghosts=True, time_restricted=plotted_interval)
