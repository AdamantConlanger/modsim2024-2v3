import matplotlib.pyplot as plt
import matplotlib.colors as clr
import math

def make_subtitle(p0, q0, x0, y0, r0, s0, k1, k2, k3, k4, is_focus):
    result = ""
    if not is_focus[0]:
        if result != "":
            result += "; "
        result += f"{p0=}"
    if not is_focus[1]:
        if result != "":
            result += "; "
        result += f"{q0=}"
    if not is_focus[2]:
        if result != "":
            result += "; "
        result += f"{x0=}"
    if not is_focus[3]:
        if result != "":
            result += "; "
        result += f"{y0=}"
    if not is_focus[4]:
        if result != "":
            result += "; "
        result += f"{r0=}"
    if not is_focus[5]:
        if result != "":
            result += "; "
        result += f"{s0=}"
    if not is_focus[6]:
        if result != "":
            result += "; "
        result += f"{k1=}"
    if not is_focus[7]:
        if result != "":
            result += "; "
        result += f"{k2=}"
    if not is_focus[8]:
        if result != "":
            result += "; "
        result += f"{k3=}"
    if not is_focus[9]:
        if result != "":
            result += "; "
        result += f"{k4=}"

    return result


def make_label(p0, q0, x0, y0, r0, s0, k1, k2, k3, k4, is_focus):
    result = ""
    if is_focus[0]:
        if result != "":
            result += "; "
        result += f"{p0=}"
    if is_focus[1]:
        if result != "":
            result += "; "
        result += f"{q0=}"
    if is_focus[2]:
        if result != "":
            result += "; "
        result += f"{x0=}"
    if is_focus[3]:
        if result != "":
            result += "; "
        result += f"{y0=}"
    if is_focus[4]:
        if result != "":
            result += "; "
        result += f"{r0=}"
    if is_focus[5]:
        if result != "":
            result += "; "
        result += f"{s0=}"
    if is_focus[6]:
        if result != "":
            result += "; "
        result += f"{k1=}"
    if is_focus[7]:
        if result != "":
            result += "; "
        result += f"{k2=}"
    if is_focus[8]:
        if result != "":
            result += "; "
        result += f"{k3=}"
    if is_focus[9]:
        if result != "":
            result += "; "
        result += f"{k4=}"

    return result


def decide_color(variables, is_focus, extrema_variables):
    focused_index = is_focus.index(True)
    focused_extrema = extrema_variables[focused_index]
    colorizing_value = (variables[focused_index] - focused_extrema[0]) / (focused_extrema[1] - focused_extrema[0])
    # blue to yellow with constant Value 1 and Saturation at 1?
    # note that blue corresponds to lower and yellow to higher values of variables
    tmp = 2/3 - colorizing_value / 2
    tmp = tmp + math.ceil(tmp) if tmp < 0 else tmp
    return clr.hsv_to_rgb([tmp, 1, 1])
    # TODO: use OKLAB or something similar and make sure the hues are evenly spaced.


def decide_color_broader(variables, is_focus, extrema_variables):
    focused_index = is_focus.index(True)
    focused_extrema = extrema_variables[focused_index]
    colorizing_value = (variables[focused_index] - focused_extrema[0]) / (focused_extrema[1] - focused_extrema[0])
    # red to yellow via blue with constant Value 1 and Saturation at 1?
    # note that red corresponds to lower and yellow to higher values of variables
    tmp = 11/12 - 3 * colorizing_value / 4
    tmp = tmp + math.ceil(tmp) if tmp < 0 else tmp
    return clr.hsv_to_rgb([tmp, 1, 1])
    # TODO: use OKLAB or something similar and make sure the hues are evenly spaced.


def decide_color_2d(variables, is_focus, extrema_variables):
    focused_index = is_focus.index(True)
    next_focused_index = is_focus[focused_index+1:].index(True)
    focused_extrema = extrema_variables[focused_index]
    next_focused_extrema = extrema_variables[next_focused_index]
    colorizing_value = (variables[focused_index] - focused_extrema[0]) / (focused_extrema[1] - focused_extrema[0])
    tmp = variables[next_focused_index] - next_focused_extrema[0]
    next_colorizing_value = tmp / (next_focused_extrema[1] - next_focused_extrema[0])
    tmp = 11/12 - 3 * colorizing_value / 4
    tmp = tmp + math.ceil(tmp) if tmp < 0 else tmp
    return clr.hsv_to_rgb([tmp, 1 - next_colorizing_value / 2, 1])
    # TODO: use OKLAB or something similar and make sure the hues are evenly spaced.
    

def visualize(f, initials_list, coefficients_list, interval_list, evaluations_list, solution_list, is_focus, extrema_variables, *, suppress_legend=False, suppress_ghosts=False, time_restricted=None, broader_colors=False, colors_2d=False):
    fig, axs = plt.subplots(3, 2, layout='constrained')
    if len(solution_list) == 0:
        the_subtitle = ""
    for n in range(len(solution_list)):
        initials_out = initials_list[n]
        coefficients_out = coefficients_list[n]
        solution_out = solution_list[n]
        t = solution_out.t
        p, q, x, y, r, s = solution_out.y
        p0, q0, x0, y0, r0, s0 = initials_out
        k1, k2, k3, k4 = coefficients_out
        if n == 1:
            the_subtitle = make_subtitle(p0, q0, x0, y0, r0, s0, k1, k2, k3, k4, is_focus)
        the_label = make_label(p0, q0, x0, y0, r0, s0, k1, k2, k3, k4, is_focus)
        if broader_colors:
            the_color = decide_color_broader([p0, q0, x0, y0, r0, s0, k1, k2, k3, k4], is_focus, extrema_variables)
        elif colors_2d:
            the_color = decide_color_2d([p0, q0, x0, y0, r0, s0, k1, k2, k3, k4], is_focus, extrema_variables)
        else:
            the_color = decide_color([p0, q0, x0, y0, r0, s0, k1, k2, k3, k4], is_focus, extrema_variables)
        axs[0, 0].plot(t, p, label=the_label, linewidth=1.5, color=the_color)
        axs[0, 1].plot(t, q, label=the_label, linewidth=1.5, color=the_color)
        axs[1, 0].plot(t, x, label=the_label, linewidth=1.5, color=the_color)
        axs[1, 1].plot(t, y, label=the_label, linewidth=1.5, color=the_color)
        axs[2, 0].plot(t, r, label=the_label, linewidth=1.5, color=the_color)
        axs[2, 1].plot(t, s, label=the_label, linewidth=1.5, color=the_color)
    if not suppress_ghosts:
        for n in range(len(solution_list)):
            initials_out = initials_list[n]
            coefficients_out = coefficients_list[n]
            solution_out = solution_list[n]
            t = solution_out.t
            p, q, x, y, r, s = solution_out.y
            p0, q0, x0, y0, r0, s0 = initials_out
            k1, k2, k3, k4 = coefficients_out
            the_label = make_label(p0, q0, x0, y0, r0, s0, k1, k2, k3, k4, is_focus)
            if broader_colors:
                the_color = decide_color_broader([p0, q0, x0, y0, r0, s0, k1, k2, k3, k4], is_focus, extrema_variables)
            elif colors_2d:
                the_color = decide_color_2d([p0, q0, x0, y0, r0, s0, k1, k2, k3, k4], is_focus, extrema_variables)
            else:
                the_color = decide_color([p0, q0, x0, y0, r0, s0, k1, k2, k3, k4], is_focus, extrema_variables)
            axs[1, 0].plot(t, y, label=the_label, alpha=0.2, linewidth=1.5, color=the_color)  # for comparison
            axs[1, 1].plot(t, x, label=the_label, alpha=0.2, linewidth=1.5, color=the_color)  # for comparison
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
    if not suppress_ghosts:
        axs[1, 0].set_ylim(min(x_min, y_min), max(x_max, y_max))
        axs[1, 1].set_ylim(min(x_min, y_min), max(x_max, y_max))
    axs[2, 0].set_ylim(r_min, r_max)
    axs[2, 1].set_ylim(s_min, s_max)
    if time_restricted is not None:
        axs[0, 0].set_xlim(time_restricted[0], time_restricted[1])
        axs[0, 1].set_xlim(time_restricted[0], time_restricted[1])
        axs[1, 0].set_xlim(time_restricted[0], time_restricted[1])
        axs[1, 1].set_xlim(time_restricted[0], time_restricted[1])
        axs[2, 0].set_xlim(time_restricted[0], time_restricted[1])
        axs[2, 1].set_xlim(time_restricted[0], time_restricted[1])
    print("graph success")
    handles, labels = axs[2, 0].get_legend_handles_labels()
    if not suppress_legend:
        fig.legend(handles, labels, mode='expand', loc='outside lower center', ncols=3, fontsize="small")
    the_title = "A graph of the full simple unreduced model."
    fig.suptitle(the_title + "\n" + the_subtitle)
    plt.show()

# TODO: make it so the names of the variables and parameters are passed along so the label makers are general.
