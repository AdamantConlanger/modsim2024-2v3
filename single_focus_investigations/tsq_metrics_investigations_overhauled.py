from cartesian_product import cartesian_product
from simulate_and_plot_metrics import simulate_and_plot_metrics
import numpy as np
import math


def f(t, u, coeffs):
    x, y = u
    mu1, mu2 = coeffs
    return np.array([
        x * x * y - x - math.sqrt(mu2 / mu1) * x + math.sqrt(mu1 * mu2),
        x - x * x * y
    ])


##################################
base_initials = [0, 0]  # list of starting values of the variables; first part of the parameters.
base_coefficients = [0, 0]  # list of coefficients for reaction speeds; second part of the parameters.
interval = (0, 1600)  # cutoff point in time to stop the simulation at, or None for the default value of 50.
granularity = 1600  # number of points in time to actually log the values at (not counting t=0),
# or None to let the solver itself decide for us.

vary_simultaneously = False  # whether to entrywise combine the variations (True) or Cartesian them (False)
multiplicative = False  # whether to apply variations multiplicatively (True) or additively (False)
variations_initials = [None, None]
# variations_coefficients = [np.linspace(0, 1, 51)[1:], np.linspace(0, 1, 51)[1:]]
variations_coefficients = [np.linspace(0, 2, 11)[1:], np.linspace(0, 2, 11)[1:]]



############################################

tmp = np.array([int(multiplicative)])
variations_initials = [[item] if isinstance(item, int | float) else item for item in variations_initials]
variations_coefficients = [[item] if isinstance(item, int | float) else item for item in variations_coefficients]
variations_initials = [tmp if item is None else np.array(item) for item in variations_initials]
variations_coefficients = [tmp if item is None else np.array(item) for item in variations_coefficients]

initials_length = len(base_initials)
coefficients_length = len(base_coefficients)
base_variables = base_initials + base_coefficients
variations_variables = variations_initials + variations_coefficients
total_length = len(variations_variables)

if vary_simultaneously:
    parallel_length = max([len(item) for item in variations_variables])
    for index in range(len(variations_variables)):
        if len(variations_variables[index]) == 1:
            variations_variables[index] = np.repeat(variations_variables[index], parallel_length)
        elif len(variations_variables[index]) != parallel_length:
            raise ValueError("Arrays aren't all of the same length or singletons")
    variations_combined = np.transpose(np.array(variations_variables))
else:
    variations_combined = cartesian_product(*variations_variables)


if multiplicative:
    separate_modified_variables = [base_variables[index] * variations_variables[index] for index in range(total_length)]
    variables_list = [np.array(base_variables) * item for item in variations_combined]
else:
    separate_modified_variables = [base_variables[index] + variations_variables[index] for index in range(total_length)]
    variables_list = [np.array(base_variables) + item for item in variations_combined]

variables_maxima = [max(item) for item in separate_modified_variables]
variables_minima = [min(item) for item in separate_modified_variables]

initials_list = [variables[:initials_length] for variables in variables_list]
coefficients_list = [variables[initials_length:] for variables in variables_list]

mu1_minimum = variables_minima[2]
mu1_maximum = variables_maxima[2]
mu2_minimum = variables_minima[3]
mu2_maximum = variables_maxima[3]
coeffs_extents = (mu1_minimum, mu1_maximum, mu2_minimum, mu2_maximum)

interval_list = [interval for item in variables_list]

if granularity is None:
    evaluations_list = None
else:
    evaluations_list = [np.linspace(interval[0], interval[1], granularity + 1) for item in variables_list]

# result = simulate(f, initials_list, coefficients_list, interval_list, evaluations_list)
# f, initials_list, coefficients_list, interval_list, evaluations_list, solution_list = result

# # freeing memory
# result = None
# f = None
# initials_list = None
# interval_list = None
# evaluations_list = None

#metrics = determine_and_plot_metrics(coefficients_list, solution_list, separate_modified_variables, coeffs_extents, (10**-4, 10**-4), use_relative=False)

simulate_and_plot_metrics(f, base_initials, coefficients_list, interval, granularity, atol=10**-6, rtol=10**-5, tolerance=(10**-4, 10**-4), evaluations_list=evaluations_list, use_relative=False,
                          plot_periods=True, plot_upper_x_amps=True, plot_lower_x_amps=True, plot_upper_y_amps=True, plot_lower_y_amps=True, plot_collapse_moments=True)
