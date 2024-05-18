def perform_program(simulate, cartesian_product):
    import numpy as np
    import matplotlib.pyplot as plt


    def make_subtitle(items, names):
        result = ""
        for index in range(len(items)):
            result += ("" if result == "" else "; ") + names[index] + f"={items[index]}"
        return result

    def visualize(item_names, initials, coefficients, solution, *, show_legend=True, show_ghosts=False, paired_bounds=True, plotted_interval=None, mini_text=False, mini_mini_text=False):
        fig, axs = plt.subplots(3, 2, layout='constrained')
        t = solution.t
        p, q, x, y, r, s = solution.y
        p0, q0 = initials[:2]
        items = list(initials) + list(coefficients)
        the_subtitle = make_subtitle(items, item_names)
        axs[0, 0].plot(t, p, linewidth=1.5)
        axs[0, 1].plot(t, q, linewidth=1.5)
        axs[1, 0].plot(t, x, linewidth=1.5)
        axs[1, 1].plot(t, y, linewidth=1.5)
        axs[2, 0].plot(t, r, linewidth=1.5)
        axs[2, 1].plot(t, s, linewidth=1.5)
        if show_ghosts or paired_bounds:
            axs[1, 0].plot(t, y, alpha=0.2, linewidth=1.5)  # for comparison
            axs[1, 1].plot(t, x, alpha=0.2, linewidth=1.5)  # for comparison
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

    item_names = ["p0", "q0", "x0", "y0", "r0", "s0", "k1", "k2", "k3", "k4"] # names of initials and coeffs.
    initials = [0.59, 0.56, 0.83, 1.25, 0, 0]  # list of starting values of the variables.
    coefficients = [0.19, 1.07, 0.85, 0.22]  # list of coefficients for reaction speeds.
    interval = (0, 500)  # cutoff point in time to stop the simulation at, or None for the default value of 50.
    granularity = 5000  # number of points in time to actually log the values at (not counting t=0),
    # or None to let the solver itself decide for us.
    plotted_interval = None  # time span to actually plot, as closed interval. or None for full plot.
    show_ghosts = False  # whether to show faint ghosts of the plots of y and x in the graphs for x and y or not.
    paired_bounds = True  # whether to force the graphs for x and y to use the same graph extent
    show_legend = False  # whether to add a legend or not.
    text_smallness = 0  # 0 for standard legend text size, 1 for smaller, 2 for tiny.
    absolute_tolerance = 10**-7  # absolute tolerance of the simulation.
    relative_tolerance = 10**-6  # relative tolerance of the simulation.

    evaluations = None if granularity is None else np.linspace(interval[0], interval[1], granularity + 1)
    mini_text = text_smallness == 1
    mini_mini_text = text_smallness == 2

    solution, = simulate(f, [initials], [coefficients], interval, evaluations, absolute_tolerance, relative_tolerance)

    visualize(item_names, initials, coefficients, solution, show_legend=show_legend, show_ghosts=show_ghosts, paired_bounds=paired_bounds,
              plotted_interval=plotted_interval, mini_text=mini_text, mini_mini_text=mini_mini_text)
