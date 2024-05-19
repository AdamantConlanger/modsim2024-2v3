import scipy.integrate as scint
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import math

plt.rcParams.update({'font.size': 18})

# https://stackoverflow.com/a/11146645/18375328


def cartesian_product(*arrays):
    broadcastable = np.ix_(*arrays)
    broadcasted = np.broadcast_arrays(*broadcastable)
    rows, cols = np.prod(broadcasted[0].shape), len(broadcasted)
    dtype = np.result_type(*arrays)

    out = np.empty(rows * cols, dtype=dtype)
    start, end = 0, rows
    for a in broadcasted:
        out[start:end] = a.reshape(-1)
        start, end = end, end + rows
    return out.reshape(cols, rows).T


def cartesian_product_for_variations(variations_initials, variations_coefficients):
    cart_prod = cartesian_product(*(variations_initials + variations_coefficients))
    len_initials = len(variations_initials)
    len_both = len_initials + len(variations_coefficients)
    return [(item[0:len_initials], item[len_initials:len_both]) for item in cart_prod]



def simulate(f, initials_list, coefficients_list, interval, evaluations=None, atol=10**-7, rtol=10**-6):
    solution_list = list()

    for iteration in range(len(initials_list)):
        sol = scint.solve_ivp(fun=f,
                              t_span=interval,
                              y0=initials_list[iteration],
                              method="Radau",
                              t_eval=evaluations,
                              args=(coefficients_list[iteration],),
                              atol=atol,
                              rtol=rtol)

        solution_list.append(sol)

    return solution_list


def make_subtitle(items, names, focus):
    result = ""
    for index in range(len(items)):
        if not focus[index]:
            result += ("" if result == "" else "; ") + names[index] + f"={items[index]}"
    return result

def make_label(items, names, focus):
    result = ""
    for index in range(len(items)):
        if focus[index]:
            result += ("" if result == "" else "; ") + names[index] + f"={items[index]}"
    return result


def decide_color_by_index(index, number_of_indices, color_stops=[], invert_colors=False):
    number_of_stops_before = [index >= item for item in color_stops].count(True)
    number_of_stops_total = len(color_stops)
    stop_size = 1/2  # size of a separating region relative to size of a color region
    total_surplus_value_from_stops = number_of_stops_total * stop_size / (number_of_stops_total + 1)
    value_from_stops = number_of_stops_before * stop_size / (number_of_stops_total + 1)
    uncorrected_colorizing_value = index / (number_of_indices - 1) if number_of_indices != 1 else 0
    colorizing_value = (uncorrected_colorizing_value + value_from_stops) / (1 + total_surplus_value_from_stops)
    print(index, number_of_stops_before, number_of_stops_total, total_surplus_value_from_stops, value_from_stops, uncorrected_colorizing_value, colorizing_value)
    if colorizing_value > 1:
        colorizing_value = 1
    if colorizing_value < 0:
        colorizing_value = 0
    if invert_colors:
        colorizing_value = 1 - colorizing_value
    # blue to yellow with constant Value 1 and Saturation at 1?
    # note that blue corresponds to lower and yellow to higher values of variables
    tmp = 2/3 - colorizing_value / 2
    tmp = tmp + math.ceil(tmp) if tmp < 0 else tmp
    return clr.hsv_to_rgb([tmp, 1, 1])
    # TODO: use OKLAB or something similar and make sure the hues are evenly spaced.


def decide_color_by_index_broader(index, number_of_indices, color_stops=[], invert_colors=False):
    number_of_stops_before = [index >= item for item in color_stops].count(True)
    number_of_stops_total = len(color_stops)
    stop_size = 1/2 # size of a separating region relative to size of a color region
    total_surplus_value_from_stops = number_of_stops_total * stop_size / (number_of_stops_total + 1)
    value_from_stops = number_of_stops_before * stop_size / (number_of_stops_total + 1)
    uncorrected_colorizing_value = index / (number_of_indices - 1) if number_of_indices != 1 else 0
    colorizing_value = (uncorrected_colorizing_value + value_from_stops) / (1 + total_surplus_value_from_stops)
    if colorizing_value > 1:
        colorizing_value = 1
    if colorizing_value < 0:
        colorizing_value = 0
    if invert_colors:
        colorizing_value = 1 - colorizing_value
    # red to yellow via blue with constant Value 1 and Saturation at 1?
    # note that red corresponds to lower and yellow to higher values of variables
    tmp = 11/12 - 3 * colorizing_value / 4
    tmp = tmp + math.ceil(tmp) if tmp < 0 else tmp
    return clr.hsv_to_rgb([tmp, 1, 1])
    # TODO: use OKLAB or something similar and make sure the hues are evenly spaced.

def decide_color(variables, focus, extrema_variables, invert_colors=False):
    focused_index = focus.index(True)
    focused_extrema = extrema_variables[focused_index]
    colorizing_value = (variables[focused_index] - focused_extrema[0]) / (focused_extrema[1] - focused_extrema[0])
    if colorizing_value > 1:
        colorizing_value = 1
    if colorizing_value < 0:
        colorizing_value = 0
    if invert_colors:
        colorizing_value = 1 - colorizing_value
    # blue to yellow with constant Value 1 and Saturation at 1?
    # note that blue corresponds to lower and yellow to higher values of variables
    tmp = 2/3 - colorizing_value / 2
    tmp = tmp + math.ceil(tmp) if tmp < 0 else tmp
    return clr.hsv_to_rgb([tmp, 1, 1])
    # TODO: use OKLAB or something similar and make sure the hues are evenly spaced.

def decide_color_broader(variables, focus, extrema_variables, invert_colors=False):
    focused_index = focus.index(True)
    focused_extrema = extrema_variables[focused_index]
    colorizing_value = (variables[focused_index] - focused_extrema[0]) / (focused_extrema[1] - focused_extrema[0])
    if colorizing_value > 1:
        colorizing_value = 1
    if colorizing_value < 0:
        colorizing_value = 0
    if invert_colors:
        colorizing_value = 1 - colorizing_value
    # red to yellow via blue with constant Value 1 and Saturation at 1?
    # note that red corresponds to lower and yellow to higher values of variables
    tmp = 11/12 - 3 * colorizing_value / 4
    tmp = tmp + math.ceil(tmp) if tmp < 0 else tmp
    return clr.hsv_to_rgb([tmp, 1, 1])
    # TODO: use OKLAB or something similar and make sure the hues are evenly spaced.

def decide_color_2d(variables, focus, extrema_variables, invert_colors=False):
    focused_index = focus.index(True)
    next_focused_index = focus[focused_index+1:].index(True) + focused_index
    focused_extrema = extrema_variables[focused_index]
    next_focused_extrema = extrema_variables[next_focused_index]
    colorizing_value = (variables[focused_index] - focused_extrema[0]) / (focused_extrema[1] - focused_extrema[0])
    if colorizing_value > 1:
        colorizing_value = 1
    if colorizing_value < 0:
        colorizing_value = 0
    if invert_colors:
        colorizing_value = 1 - colorizing_value
    tmp = variables[next_focused_index] - next_focused_extrema[0]
    next_colorizing_value = tmp / (next_focused_extrema[1] - next_focused_extrema[0])
    if next_colorizing_value > 1:
        next_colorizing_value = 1
    if next_colorizing_value < 0:
        next_colorizing_value = 0
    if invert_colors:
        next_colorizing_value = 1 - next_colorizing_value
    tmp = 11/12 - 3 * colorizing_value / 4
    tmp = tmp + math.ceil(tmp) if tmp < 0 else tmp
    tmp2 = 1 - 5 * (1 - colorizing_value**2) * next_colorizing_value / 6
    tmp3 = 1 - colorizing_value**2 * next_colorizing_value / 4
    return clr.hsv_to_rgb([tmp, tmp2, tmp3])
    # TODO: use OKLAB or something similar and make sure the hues are evenly spaced.

def visualize(item_names, initials_list, coefficients_list, solution_list, focus, extrema_variables, *, show_legend=True, show_ghosts=False, paired_bounds=True, plotted_interval=None, broader_colors=False, colors_2d=False, mini_text=False, mini_mini_text=False, linewidth=1.5, invert_colors=False, color_stops=[]):
    fig, axs = plt.subplots(1, 2, layout='constrained')
    if len(solution_list) == 0:
        the_subtitle = ""
    for n in range(len(solution_list)):
        initials_out = initials_list[n]
        coefficients_out = coefficients_list[n]
        solution_out = solution_list[n]
        t = solution_out.t
        x, y = solution_out.y
        items_out = list(initials_out) + list(coefficients_out)
        if n == 1:
            the_subtitle = make_subtitle(items_out, item_names, focus)
        the_label = make_label(items_out, item_names, focus)
        if broader_colors:
            the_color = decide_color_by_index_broader(
                n, len(solution_list), color_stops=color_stops, invert_colors=invert_colors)
        elif colors_2d:
            the_color = decide_color_2d(items_out, focus, extrema_variables, invert_colors=invert_colors)
        else:
            print(color_stops)
            the_color = decide_color_by_index(
                n, len(solution_list), color_stops=color_stops, invert_colors=invert_colors)
        axs[0].plot(t, x, label=the_label, linewidth=linewidth, color=the_color)
        axs[1].plot(t, y, label=the_label, linewidth=linewidth, color=the_color)
        if show_ghosts:
            axs[0].plot(t, y, label=the_label, alpha=0.2, linewidth=linewidth, color=the_color)  # for comparison
            axs[1].plot(t, x, label=the_label, alpha=0.2, linewidth=linewidth, color=the_color)  # for comparison
            # TODO: make the color of these lines appear in the legend too
            # TODO: make these lines appear behind the other ones, but with these colors
    # TODO: make it so the labels are aligned with one another
    axs[0].set(ylabel="x", xlabel="t")
    axs[0].grid(True, linestyle='dashed')
    axs[1].set(ylabel="y", xlabel="t")
    axs[1].grid(True, linestyle='dashed')
    x_min, x_max = axs[0].get_ylim()
    y_min, y_max = axs[1].get_ylim()
    if show_ghosts or paired_bounds:
        axs[0].set_ylim(min(x_min, y_min), max(x_max, y_max))
        axs[1].set_ylim(min(x_min, y_min), max(x_max, y_max))
    if plotted_interval is not None:
        axs[0].set_xlim(plotted_interval[0], plotted_interval[1])
        axs[1].set_xlim(plotted_interval[0], plotted_interval[1])
    print("graph success")
    handles, labels = axs[0].get_legend_handles_labels()
    if show_legend:
        fontsize = "xx-small" if mini_mini_text else "x-small" if mini_text else "small"
        fig.legend(handles, labels, mode='expand', loc='outside lower center', ncols=5, fontsize=fontsize)
    the_title = "A graph of the simple truncated unreduced model."
    fig.suptitle(the_title + "\n" + the_subtitle)
    plt.show()

# TODO: make it so the names of the variables and parameters are passed along so the label makers are general.

def f(t, u, coeffs):
    x, y = u
    k1, k2, k3, k4, p, q = coeffs
    return np.array([
        k1 * p - k2 * q * x + k3 * x * x * y - k4 * x,
        k2 * q * x - k3 * x * x * y
    ])

item_names = ["x0", "y0", "k1", "k2", "k3", "k4", "p", "q"]  # names of initials and coeffs.
base_initials = [0, 0]  # list of starting values of the variables.
base_coefficients = [0.19, 1.07, 0.85, 0, 0.59, 0.56]  # list of coefficients for reaction speeds.
interval = (0, 250)  # cutoff point in time to stop the simulation at, or None for the default value of 50.
granularity = 5000  # number of points in time to actually log the values at (not counting t=0),
# or None to let the solver itself decide for us.
plotted_interval = None  # time span to actually plot, as closed interval. or None for full plot.
show_ghosts = False  # whether to show faint ghosts of the plots of y and x in the graphs for x and y or not.
paired_bounds = False  # whether to force the graphs for x and y to use the same graph extent
show_legend = True  # whether to add a legend or not.
broader_colors = True  # whether to use a larger-than-usual color spectrum.
text_smallness = 0  # 0 for standard legend text size, 1 for smaller, 2 for tiny.
linewidth = 2  # width of plotted lines
invert_colors = False  # whether to use invert the color scheme. "False" uses blue for low values.
absolute_tolerance = 10**-7  # absolute tolerance of the simulation.
relative_tolerance = 10**-6  # relative tolerance of the simulation.
vary_simultaneously = False  # whether to entrywise combine the variations (True) or Cartesian them (False).
multiplicative = False  # whether to apply variations multiplicatively (True) or additively (False).
variations_initials = [None, None]  # variations in the initials.
# variations in the coeffs.
# my_tmp = np.concatenate((np.linspace(50, 250, 9), np.array([1000]))) / 100
# np.linspace(0, 20, 11)[1:]
# np.linspace(48, 60, 13)
# np.linspace(60, 80, 6)
my_tmp = np.concatenate((np.linspace(0, 10, 5)[1:], np.linspace(10, 22, 4)[1:], np.linspace(50, 70, 11), np.array([80]))) / 100
# my_tmp = np.concatenate((np.linspace(5, 25, 9), np.array([100]))) / 10
variations_coefficients = [None, None, None, my_tmp, None, None]
focus_initials = [False, False]  # which variations should determine plot colors?
focus_coefficients = [False, False, False, True, False, False]  # which variations should determine plot colors?

initials_length = len(base_initials)
base_variables = base_initials + base_coefficients
variations_variables = variations_initials + variations_coefficients
focus = focus_initials + focus_coefficients

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

# extrema_variables = [(0.0, 0.0), (0.0, 0.0), (0.19, 0.19), (1.07, 1.07),
#                      (0.5, 2.5), (0.22, 0.22), (0.59, 0.59), (0.56, 0.56)]

color_stops = (7,)

initials_list = [variables[:initials_length] for variables in variables_list]
coefficients_list = [variables[initials_length:] for variables in variables_list]
evaluations = None if granularity is None else np.linspace(interval[0], interval[1], granularity + 1)
colors_2d = (focus_initials + focus_coefficients).count(True) == 2
mini_text = text_smallness == 1
mini_mini_text = text_smallness == 2

solution_list = simulate(f, initials_list, coefficients_list, interval,
                            evaluations, absolute_tolerance, relative_tolerance)

visualize(item_names, initials_list, coefficients_list, solution_list, focus, extrema_variables, show_legend=show_legend,
            show_ghosts=show_ghosts, paired_bounds=paired_bounds, plotted_interval=plotted_interval, broader_colors=broader_colors, colors_2d=colors_2d, mini_text=mini_text, mini_mini_text=mini_mini_text, linewidth=linewidth, invert_colors=invert_colors, color_stops=color_stops)
