import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scint

plt.rcParams.update({'font.size': 18})

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

def make_subtitle(items, names):
    result = ""
    for index in range(len(items)):
        result += ("" if result == "" else "; ") + names[index] + f"={items[index]}"
    return result

def visualize(item_names, initials, coefficients, solution, *, show_legend=True, show_ghosts=False, paired_bounds=True, plotted_interval=None, mini_text=False, mini_mini_text=False, linewidth=1.5):
    fig, axs = plt.subplots(1, 2, layout='constrained')
    t = solution.t
    x, y = solution.y
    items = list(initials) + list(coefficients)
    the_subtitle = make_subtitle(items, item_names)
    axs[0].plot(t, x, linewidth=linewidth)
    axs[1].plot(t, y, linewidth=linewidth)
    if show_ghosts:
        axs[0].plot(t, y, alpha=0.2, linewidth=linewidth)  # for comparison
        axs[1].plot(t, x, alpha=0.2, linewidth=linewidth)  # for comparison
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

def f(t, u, coeffs):
    x, y = u
    k1, k2, k3, k4, p, q = coeffs
    return np.array([
        k1 * p - k2 * q * x + k3 * x * x * y - k4 * x,
        k2 * q * x - k3 * x * x * y
    ])

item_names = ["x0", "y0", "k1", "k2", "k3", "k4", "p", "q"]  # names of initials and coeffs.
initials = [0.83, 1.25]  # list of starting values of the variables.
coefficients = [0.19, 1.07, 0.85, 0.22, 0.59, 0.56]  # list of coefficients for reaction speeds.
interval = (0, 200)  # cutoff point in time to stop the simulation at, or None for the default value of 50.
granularity = 5000  # number of points in time to actually log the values at (not counting t=0),
# or None to let the solver itself decide for us.
plotted_interval = None  # time span to actually plot, as closed interval. or None for full plot.
show_ghosts = False  # whether to show faint ghosts of the plots of y and x in the graphs for x and y or not.
paired_bounds = True  # whether to force the graphs for x and y to use the same graph extent
show_legend = False  # whether to add a legend or not.
text_smallness = 0  # 0 for standard legend text size, 1 for smaller, 2 for tiny.
linewidth = 1.5  # width of plotted lines
absolute_tolerance = 10**-7  # absolute tolerance of the simulation.
relative_tolerance = 10**-6  # relative tolerance of the simulation.

evaluations = None if granularity is None else np.linspace(interval[0], interval[1], granularity + 1)
mini_text = text_smallness == 1
mini_mini_text = text_smallness == 2

solution, = simulate(f, [initials], [coefficients], interval, evaluations, absolute_tolerance, relative_tolerance)

visualize(item_names, initials, coefficients, solution, show_legend=show_legend, show_ghosts=show_ghosts, paired_bounds=paired_bounds,
            plotted_interval=plotted_interval, mini_text=mini_text, mini_mini_text=mini_mini_text, linewidth=linewidth)
