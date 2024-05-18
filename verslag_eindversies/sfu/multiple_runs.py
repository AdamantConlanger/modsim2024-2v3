def perform_program(simulate, cartesian_product):
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as clr
    import math

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

    def decide_color(variables, focus, extrema_variables):
        focused_index = focus.index(True)
        focused_extrema = extrema_variables[focused_index]
        colorizing_value = (variables[focused_index] - focused_extrema[0]) / (focused_extrema[1] - focused_extrema[0])
        # blue to yellow with constant Value 1 and Saturation at 1?
        # note that blue corresponds to lower and yellow to higher values of variables
        tmp = 2/3 - colorizing_value / 2
        tmp = tmp + math.ceil(tmp) if tmp < 0 else tmp
        return clr.hsv_to_rgb([tmp, 1, 1])
        # TODO: use OKLAB or something similar and make sure the hues are evenly spaced.

    def decide_color_broader(variables, focus, extrema_variables):
        focused_index = focus.index(True)
        focused_extrema = extrema_variables[focused_index]
        colorizing_value = (variables[focused_index] - focused_extrema[0]) / (focused_extrema[1] - focused_extrema[0])
        # red to yellow via blue with constant Value 1 and Saturation at 1?
        # note that red corresponds to lower and yellow to higher values of variables
        tmp = 11/12 - 3 * colorizing_value / 4
        tmp = tmp + math.ceil(tmp) if tmp < 0 else tmp
        return clr.hsv_to_rgb([tmp, 1, 1])
        # TODO: use OKLAB or something similar and make sure the hues are evenly spaced.

    def decide_color_2d(variables, focus, extrema_variables):
        focused_index = focus.index(True)
        next_focused_index = focus[focused_index+1:].index(True) + focused_index
        focused_extrema = extrema_variables[focused_index]
        next_focused_extrema = extrema_variables[next_focused_index]
        colorizing_value = (variables[focused_index] - focused_extrema[0]) / (focused_extrema[1] - focused_extrema[0])
        tmp = variables[next_focused_index] - next_focused_extrema[0]
        next_colorizing_value = tmp / (next_focused_extrema[1] - next_focused_extrema[0])
        tmp = 11/12 - 3 * colorizing_value / 4
        tmp = tmp + math.ceil(tmp) if tmp < 0 else tmp
        tmp2 = 1 - 5 * (1 - colorizing_value**2) * next_colorizing_value / 6
        tmp3 = 1 - colorizing_value**2 * next_colorizing_value / 4
        return clr.hsv_to_rgb([tmp, tmp2, tmp3])
        # TODO: use OKLAB or something similar and make sure the hues are evenly spaced.

    def visualize(item_names, initials_list, coefficients_list, solution_list, focus, extrema_variables, *, show_legend=True, show_ghosts=False, paired_bounds=True, plotted_interval=None, broader_colors=False, colors_2d=False, mini_text=False, mini_mini_text=False, linewidth=1.5):
        fig, axs = plt.subplots(3, 2, layout='constrained')
        if len(solution_list) == 0:
            the_subtitle = ""
        for n in range(len(solution_list)):
            initials_out = initials_list[n]
            coefficients_out = coefficients_list[n]
            solution_out = solution_list[n]
            t = solution_out.t
            p, q, x, y, r, s = solution_out.y
            p0, q0 = initials_out[:2]
            items_out = list(initials_out) + list(coefficients_out)
            if n == 1:
                the_subtitle = make_subtitle(items_out, item_names, focus)
            the_label = make_label(items_out, item_names, focus)
            if broader_colors:
                the_color = decide_color_broader(items_out, focus, extrema_variables)
            elif colors_2d:
                the_color = decide_color_2d(items_out, focus, extrema_variables)
            else:
                the_color = decide_color(items_out, focus, extrema_variables)
            axs[0, 0].plot(t, p, label=the_label, linewidth=linewidth, color=the_color)
            axs[0, 1].plot(t, q, label=the_label, linewidth=linewidth, color=the_color)
            axs[1, 0].plot(t, x, label=the_label, linewidth=linewidth, color=the_color)
            axs[1, 1].plot(t, y, label=the_label, linewidth=linewidth, color=the_color)
            axs[2, 0].plot(t, r, label=the_label, linewidth=linewidth, color=the_color)
            axs[2, 1].plot(t, s, label=the_label, linewidth=linewidth, color=the_color)
            if show_ghosts:
                axs[1, 0].plot(t, y, label=the_label, alpha=0.2, linewidth=linewidth, color=the_color)  # for comparison
                axs[1, 1].plot(t, x, label=the_label, alpha=0.2, linewidth=linewidth, color=the_color)  # for comparison
                # TODO: make the color of these lines appear in the legend too
                # TODO: make these lines appear behind the other ones, but with these colors
        # TODO: make it so the labels are aligned with one another
        axs[0, 0].set(ylabel="p")
        axs[0, 0].xaxis.set_major_formatter('')
        axs[0, 0].tick_params(axis='x', length=0)
        axs[0, 0].grid(True, linestyle='dashed')
        axs[0, 1].set(ylabel="q")
        axs[0, 1].xaxis.set_major_formatter('')
        axs[0, 1].tick_params(axis='x', length=0)
        axs[0, 1].grid(True, linestyle='dashed')
        axs[1, 0].set(ylabel="x")
        axs[1, 0].xaxis.set_major_formatter('')
        axs[1, 0].tick_params(axis='x', length=0)
        axs[1, 0].grid(True, linestyle='dashed')
        axs[1, 1].set(ylabel="y")
        axs[1, 1].xaxis.set_major_formatter('')
        axs[1, 1].tick_params(axis='x', length=0)
        axs[1, 1].grid(True, linestyle='dashed')
        axs[2, 0].set(ylabel="r", xlabel='t')
        axs[2, 0].grid(True, linestyle='dashed')
        axs[2, 1].set(ylabel="s", xlabel='t')
        axs[2, 1].grid(True, linestyle='dashed')
        p_min, p_max = axs[0, 0].get_ylim()
        q_min, q_max = axs[0, 1].get_ylim()
        x_min, x_max = axs[1, 0].get_ylim()
        y_min, y_max = axs[1, 1].get_ylim()
        r_min, r_max = axs[2, 0].get_ylim()
        s_min, s_max = axs[2, 1].get_ylim()
        axs[0, 0].set_ylim(0, 2 * p0)
        axs[0, 1].set_ylim(0, 2 * q0)
        if show_ghosts or paired_bounds:
            axs[1, 0].set_ylim(min(x_min, y_min), max(x_max, y_max))
            axs[1, 1].set_ylim(min(x_min, y_min), max(x_max, y_max))
        axs[2, 0].set_ylim(r_min, r_max)
        axs[2, 1].set_ylim(s_min, s_max)
        if plotted_interval is not None:
            axs[0, 0].set_xlim(plotted_interval[0], plotted_interval[1])
            axs[0, 1].set_xlim(plotted_interval[0], plotted_interval[1])
            axs[1, 0].set_xlim(plotted_interval[0], plotted_interval[1])
            axs[1, 1].set_xlim(plotted_interval[0], plotted_interval[1])
            axs[2, 0].set_xlim(plotted_interval[0], plotted_interval[1])
            axs[2, 1].set_xlim(plotted_interval[0], plotted_interval[1])
        print("graph success")
        handles, labels = axs[2, 0].get_legend_handles_labels()
        if show_legend:
            fontsize = "xx-small" if mini_mini_text else "x-small" if mini_text else "small"
            fig.legend(handles, labels, mode='expand', loc='outside lower center', ncols=5, fontsize=fontsize)
        the_title = "A graph of the simple full unreduced model."
        fig.suptitle(the_title + "\n" + the_subtitle)
        plt.show()

    # TODO: make it so the names of the variables and parameters are passed along so the label makers are general.

    def f(t, u, coeffs):
        p, q, x, y, r, s = u
        k1, k2, k3, k4 = coeffs
        return np.array([
            0,
            0,
            k1 * p - k2 * q * x + k3 * x * x * x * y - k4 * x,
            k2 * q * x - k3 * x * x * x * y,
            k2 * q * x,
            k4 * x
        ])

    item_names = ["p0", "q0", "x0", "y0", "r0", "s0", "k1", "k2", "k3", "k4"]  # names of initials and coeffs.
    base_initials = [0.59, 0.56, 0.83, 1.25, 0, 0]  # list of starting values of the variables.
    base_coefficients = [0.19, 1.07, 0.85, 0.22]  # list of coefficients for reaction speeds.
    interval = (0, 1000)  # cutoff point in time to stop the simulation at, or None for the default value of 50.
    granularity = 5000  # number of points in time to actually log the values at (not counting t=0),
    # or None to let the solver itself decide for us.
    plotted_interval = None  # time span to actually plot, as closed interval. or None for full plot.
    show_ghosts = False  # whether to show faint ghosts of the plots of y and x in the graphs for x and y or not.
    paired_bounds = True  # whether to force the graphs for x and y to use the same graph extent
    show_legend = True  # whether to add a legend or not.
    broader_colors = False  # whether to use a larger-than-usual color spectrum.
    text_smallness = 0  # 0 for standard legend text size, 1 for smaller, 2 for tiny.
    linewidth = 1.5  # width of plotted lines
    absolute_tolerance = 10**-7  # absolute tolerance of the simulation.
    relative_tolerance = 10**-6  # relative tolerance of the simulation.
    vary_simultaneously = False  # whether to entrywise combine the variations (True) or Cartesian them (False).
    multiplicative = True  # whether to apply variations multiplicatively (True) or additively (False).
    variations_initials = [None, None, None, None, None, None]  # variations in the initials.
    variations_coefficients = [np.linspace(1.3, 1.35, 20), None, None, None]  # variations in the coeffs.
    focus_initials = [False, False, False, False, False, False]  # which variations should determine plot colors?
    focus_coefficients = [True, False, False, False]  # which variations should determine plot colors?

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

    initials_list = [variables[:initials_length] for variables in variables_list]
    coefficients_list = [variables[initials_length:] for variables in variables_list]
    evaluations = None if granularity is None else np.linspace(interval[0], interval[1], granularity + 1)
    colors_2d = (focus_initials + focus_coefficients).count(True) == 2
    mini_text = text_smallness == 1
    mini_mini_text = text_smallness == 2

    solution_list = simulate(f, initials_list, coefficients_list, interval,
                             evaluations, absolute_tolerance, relative_tolerance)

    visualize(item_names, initials_list, coefficients_list, solution_list, focus, extrema_variables, show_legend=show_legend,
              show_ghosts=show_ghosts, paired_bounds=paired_bounds, plotted_interval=plotted_interval, broader_colors=broader_colors, colors_2d=colors_2d, mini_text=mini_text, mini_mini_text=mini_mini_text, linewidth=linewidth)
