import matplotlib.pyplot as plt

def visualize_simple_unreduced_2(f, initials_list, coefficients_list, interval_list, evaluations_list, solution_list, *, suppress_legend=False, suppress_ghosts=False):
    fig, axs = plt.subplots(1, 2, layout='constrained')
    for n in range(len(solution_list)):
        initials_out = initials_list[n]
        coefficients_out = coefficients_list[n]
        solution_out = solution_list[n]
        t = solution_out.t
        p, q, x, y, r, s = solution_out.y
        p0, q0, x0, y0, r0, s0 = initials_out
        k1, k2, k3, k4 = coefficients_out
        # TODO: make it so only the varied parameters are shown here, and make the legend choose orientation
        the_label = f"u0=({p0}, {q0}, {x0}, {y0}, {r0}, {s0});\n{k1=}; {k2=}; {k3=}; {k4=}"
        axs[0].plot(t, x, label=the_label, linewidth=1.5)
        axs[1].plot(t, y, label=the_label, linewidth=1.5)
    if not suppress_ghosts:
        for n in range(len(solution_list)):
            initials_out = initials_list[n]
            coefficients_out = coefficients_list[n]
            solution_out = solution_list[n]
            t = solution_out.t
            p, q, x, y, r, s = solution_out.y
            p0, q0, x0, y0, r0, s0 = initials_out
            k1, k2, k3, k4 = coefficients_out
            the_label = f"u0=({p0}, {q0}, {x0}, {y0}, {r0}, {s0});\n{k1=}; {k2=}; {k3=}; {k4=}"
            axs[0].plot(t, y, label=the_label, alpha=0.2, linewidth=1.5)  # for comparison
            axs[1].plot(t, x, label=the_label, alpha=0.2, linewidth=1.5)  # for comparison
            # TODO: make the color of these lines appear in the legend too
            # TODO: make these lines appear behind the other ones, but with these colors
    # TODO: make it so the labels are aligned with one another
    axs[0].set(ylabel="x", xlabel='t')
    axs[1].set(ylabel="y", xlabel='t')
    axs[1].grid(True, linestyle='dashed')
    x_min, x_max = axs[0].get_ylim()
    y_min, y_max = axs[1].get_ylim()
    axs[0].set_ylim(min(x_min, y_min), max(x_max, y_max))
    axs[1].set_ylim(min(x_min, y_min), max(x_max, y_max))
    print("graph success")
    handles, labels = axs[0].get_legend_handles_labels()
    if not suppress_legend:
        fig.legend(handles, labels, mode='expand', loc='outside lower center')
    fig.suptitle("A graph of the unreduced base model with constant P and Q.")
    plt.show()
