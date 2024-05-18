def perform_program(simulate, cartesian_product):
    import numpy as np
    import matplotlib.pyplot as plt
    import math


    def make_subtitle(items, names):
        result = ""
        for index in range(len(items)):
            result += ("" if result == "" else "; ") + names[index] + f"={items[index]}"
        return result

    def visualize(item_names, initials, coefficients, solution, *, show_legend=True, show_ghosts=False, paired_bounds=True, plotted_interval=None, mini_text=False, mini_mini_text=False):
        fig, axs = plt.subplots(1, 2, layout='constrained')
        t = solution.t
        x, y = solution.y
        items = list(initials) + list(coefficients)
        the_subtitle = make_subtitle(items, item_names)
        axs[0].plot(t, x, linewidth=1.5)
        axs[1].plot(t, y, linewidth=1.5)
        if show_ghosts:
            axs[0].plot(t, y, alpha=0.2, linewidth=1.5)  # for comparison
            axs[1].plot(t, x, alpha=0.2, linewidth=1.5)  # for comparison
            # TODO: make the color of these lines appear in the legend too
            # TODO: make these lines appear behind the other ones, but with these colors
        # TODO: make it so the labels are aligned with one another
        axs[0].set(ylabel="reduced x", xlabel="reduced t")
        axs[0].grid(True, linestyle='dashed')
        axs[1].set(ylabel="reduced y", xlabel="reduced t")
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
        the_title = "A graph of the cubic truncated reduced model."
        fig.suptitle(the_title + "\n" + the_subtitle)
        plt.show()

    # TODO: make it so the names of the variables and parameters are passed along so the label makers are general.

    def f(t, u, coeffs):
        x, y = u
        mu1, mu2 = coeffs
        return np.array([
            x * x * x * y - x - math.sqrt(mu2 / mu1) * x + math.sqrt(mu2 * mu1),
            x - x * x * x * y
        ])

    item_names = ["reduced x0", "reduced y0", "alpha/beta", "alpha*beta"] # names of initials and coeffs.
    initials = [0, 0]  # list of starting values of the variables.
    coefficients = [2/3, 3/50]  # list of coefficients for reaction speeds.
    interval = (0, 3200)  # cutoff point in time to stop the simulation at, or None for the default value of 50.
    granularity = 1600  # number of points in time to actually log the values at (not counting t=0),
    # or None to let the solver itself decide for us.
    plotted_interval = None  # time span to actually plot, as closed interval. or None for full plot.
    show_ghosts = False  # whether to show faint ghosts of the plots of y and x in the graphs for x and y or not.
    paired_bounds = True  # whether to force the graphs for x and y to use the same graph extent
    show_legend = False # whether to add a legend or not.
    text_smallness = 0 # 0 for standard legend text size, 1 for smaller, 2 for tiny.
    absolute_tolerance = 10**-7 # absolute tolerance of the simulation.
    relative_tolerance = 10**-6 # relative tolerance of the simulation.

    evaluations = None if granularity is None else np.linspace(interval[0], interval[1], granularity + 1)
    mini_text = text_smallness == 1
    mini_mini_text = text_smallness == 2

    solution, = simulate(f, [initials], [coefficients], interval, evaluations, absolute_tolerance, relative_tolerance)

    visualize(item_names, initials, coefficients, solution, show_legend=show_legend, show_ghosts=show_ghosts, paired_bounds=paired_bounds, plotted_interval=plotted_interval, mini_text=mini_text, mini_mini_text=mini_mini_text)
