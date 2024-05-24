from simulate import simulate
from cartesian_product import cartesian_product
from determine_period import determine_metrics
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import numpy as np
import math


def f(t, u, coeffs):
    x, y = u
    mu1, mu2 = coeffs
    return np.array([
        x * x * y - x - math.sqrt(mu2 / mu1) * x + math.sqrt(mu1 * mu2),
        x - x * x * y
    ])


base_initials = [1, 1]  # list of starting values of the variables; first part of the parameters.
base_coefficients = [0.1, 0]  # list of coefficients for reaction speeds; second part of the parameters.
interval = (0, 400)  # cutoff point in time to stop the simulation at, or None for the default value of 50.
granularity = 1600  # number of points in time to actually log the values at (not counting t=0),
# or None to let the solver itself decide for us.
plotted_interval = None  # time span to actually plot, as closed interval. or None for full plot.

vary_simultaneously = False  # whether to entrywise combine the variations (True) or Cartesian them (False)
multiplicative = False  # whether to apply variations multiplicatively (True) or additively (False)
variations_initials = [None, None]
variations_coefficients = [np.linspace(0, 0.2, 11)[1:], np.linspace(0, 0.1, 11)[1:]]
is_focus_initials = [False, False]
is_focus_coefficients = [True, True]

initials_length = len(base_initials)
coefficients_length = len(base_coefficients)
base_variables = base_initials + base_coefficients
variations_variables = variations_initials + variations_coefficients
is_focus = is_focus_initials + is_focus_coefficients

tmp = np.array([int(multiplicative)])
variations_variables = [[item] if isinstance(item, int | float) else item for item in variations_variables]
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

# freeing memory
result = None

#visualize(f, initials_list, coefficients_list, interval_list,
#          evaluations_list, solution_list, is_focus, extrema_variables, suppress_legend=False, suppress_ghosts=True, time_restricted=plotted_interval, broader_colors=True, colors_2d=False)

metrics = determine_metrics(f, initials_list, coefficients_list, interval_list, evaluations_list, solution_list, (10**-4, 10**-4), use_relative=False)

# freeing memory
f = None
initials_list = None
interval_list = None
evaluations_list = None
solution_list = None

# coefficients_to_indices = dict()
first_coeffs = variations_coefficients[0]
second_coeffs = variations_coefficients[1]
len_first = len(first_coeffs)
len_second = len(second_coeffs)
# for first_index in range(len_first):
#     first_coeff = first_coeffs[first_index]
#     for second_index in range(len_second):
#         second_coeff = second_coeffs[second_index]
#         coefficients_to_indices[(first_coeff, second_coeff)] = first_index * len_second + second_index

periods = [item["estimated_period"] if "estimated_period" in item else 0 for item in metrics]

# freeing memory
metrics = None

# gridded_data = [[periods[first_index * len_second + second_index]
#                  for second_index in range(len_second)] for first_index in range(len_first)]
# spaced_data = [(first_coeffs[first_index], second_coeffs[second_index],
#                    periods[first_index * len_second + second_index]) for second_index in range(len_second) for first_index in range(len_first)]

# spaced_data_x = [item[0] for item in spaced_data]
# spaced_data_y = [item[1] for item in spaced_data]
# spaced_data_z = [item[2] for item in spaced_data]

first_max = max(first_coeffs)
first_min = min(first_coeffs)
second_max = max(second_coeffs)
second_min = min(second_coeffs)
desired_range_x, desired_range_y = np.meshgrid(np.linspace(first_min, first_max, len_first),
                                               np.linspace(second_min, second_max, len_second), indexing='ij')
# print(f"{first_coeffs=}")
# print(f"{second_coeffs=}")
# print(f"{desired_range_x=}")
# print(f"{desired_range_y=}")
# print(f"{coefficients_list=}")
desired_range = (desired_range_x, desired_range_y)
current_points = coefficients_list

# https://stackoverflow.com/a/14140554/18375328,
# which is superseded by https://docs.scipy.org/doc/scipy/tutorial/interpolate/ND_unstructured.html

plottable_data = griddata(coefficients_list, periods, desired_range, method="nearest")


fig, axs = plt.subplots(1, 1, layout='constrained')

use_smooth_block_interpolation = False
use_smooth_lossy_interpolation = False

extent = (first_min, first_max, second_min, second_max)

if use_smooth_block_interpolation:
    im = axs.imshow(plottable_data.T, extent=extent, interpolation="spline16", origin='lower', cmap="viridis")
elif use_smooth_lossy_interpolation:
    im = axs.imshow(plottable_data.T, extent=extent, interpolation="bicubic", origin='lower', cmap="viridis")
else:
    im = axs.imshow(plottable_data.T, extent=extent, interpolation="none", origin='lower', cmap="viridis")
axs.set_ylabel("mu2")
axs.set_xlabel("mu1")
fig.suptitle("A graph of the period of oscillation for the truncated simple (reduced) quoduct model.\nx0=1, y0=1.")
fig.colorbar(im, ax=axs, label='period of oscillation (0 if n/a)')
plt.show()

# TODO: merge simulation and metric-determining parts so we don't need 4GB or something to simulate it.
# TODO: make it so uncertain values are left out