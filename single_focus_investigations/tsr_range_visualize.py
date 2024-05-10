import matplotlib.pyplot as plt
import matplotlib.colors as clr

def make_subtitle(x0, y0, kappa1, kappa2, is_focus):
    result = ""
    if not is_focus[0]:
        if result != "":
            result += "; "
        result += f"{x0=}"
    if not is_focus[1]:
        if result != "":
            result += "; "
        result += f"{y0=}"
    if not is_focus[2]:
        if result != "":
            result += "; "
        result += f"{kappa1=}"
    if not is_focus[3]:
        if result != "":
            result += "; "
        result += f"{kappa2=}"

    return result


def make_label(x0, y0, kappa1, kappa2, is_focus):
    result = ""
    if is_focus[0]:
        if result != "":
            result += "; "
        result += f"{x0=}"
    if is_focus[1]:
        if result != "":
            result += "; "
        result += f"{y0=}"
    if is_focus[2]:
        if result != "":
            result += "; "
        result += f"{kappa1=}"
    if is_focus[3]:
        if result != "":
            result += "; "
        result += f"{kappa2=}"

    return result


def decide_color(variables, is_focus, extrema_variables):
    focused_index = is_focus.index(True)
    focused_extrema = extrema_variables[focused_index]
    colorizing_value = (variables[focused_index] - focused_extrema[0]) / (focused_extrema[1] - focused_extrema[0])
    # blue to yellow with constant Value 1 and Saturation at 1?
    # note that blue corresponds to lower and yellow to higher values
    return clr.hsv_to_rgb([2/3 - colorizing_value / 2, 1, 1])
    # TODO: use OKLAB or something similar and make sure the hues are evenly spaced.
    

def visualize(f, initials_list, coefficients_list, interval_list, evaluations_list, solution_list, is_focus, extrema_variables, *, suppress_legend=False, suppress_ghosts=False):
    fig, axs = plt.subplots(1, 2, layout='constrained')
    if len(solution_list) == 0:
        the_subtitle = ""
    for n in range(len(solution_list)):
        initials_out = initials_list[n]
        coefficients_out = coefficients_list[n]
        solution_out = solution_list[n]
        t = solution_out.t
        x, y = solution_out.y
        x0, y0 = initials_out
        kappa1, kappa2 = coefficients_out
        if n == 1:
            the_subtitle = make_subtitle(x0, y0, kappa1, kappa2, is_focus)
        the_label = make_label(x0, y0, kappa1, kappa2, is_focus)
        the_color = decide_color([x0, y0, kappa1, kappa2], is_focus, extrema_variables)
        axs[0].plot(t, x, label=the_label, linewidth=1.5, color=the_color)
        axs[1].plot(t, y, label=the_label, linewidth=1.5, color=the_color)
    if not suppress_ghosts:
        for n in range(len(solution_list)):
            initials_out = initials_list[n]
            coefficients_out = coefficients_list[n]
            solution_out = solution_list[n]
            t = solution_out.t
            x, y = solution_out.y
            x0, y0 = initials_out
            kappa1, kappa2 = coefficients_out
            the_label = make_label(x0, y0, kappa1, kappa2, is_focus)
            the_color = decide_color([x0, y0, kappa1, kappa2], is_focus, extrema_variables)
            axs[0].plot(t, y, label=the_label, alpha=0.2, linewidth=1.5, color=the_color)  # for comparison
            axs[1].plot(t, x, label=the_label, alpha=0.2, linewidth=1.5, color=the_color)  # for comparison
            # TODO: make the color of these lines appear in the legend too
            # TODO: make these lines appear behind the other ones, but with these colors
    # TODO: make it so the labels are aligned with one another
    axs[0].set(ylabel="reduced x", xlabel='reduced t')
    axs[0].grid(True, linestyle='dashed')
    axs[1].set(ylabel="reduced y", xlabel='reduced t')
    axs[1].grid(True, linestyle='dashed')
    x_min, x_max = axs[0].get_ylim()
    y_min, y_max = axs[1].get_ylim()
    if not suppress_ghosts:
        axs[0].set_ylim(min(x_min, y_min), max(x_max, y_max))
        axs[1].set_ylim(min(x_min, y_min), max(x_max, y_max))
    print("graph success")
    handles, labels = axs[0].get_legend_handles_labels()
    if not suppress_legend:
        fig.legend(handles, labels, mode='expand', loc='outside lower center', ncols=3, fontsize="small")
    the_title = "A graph of the truncated simple reduced model."
    fig.suptitle(the_title + "\n" + the_subtitle)
    plt.show()


# TODO: make it so the names of the variables and parameters are passed along so the label makers are general.