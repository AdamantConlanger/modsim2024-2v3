from new_simulator import simulate
import numpy as np
import matplotlib.pyplot as plt


def f(t, u, coeffs):
    p, q, x, y, r, s = u
    k1, k2, k3, k4 = coeffs
    return np.array([
        0,
        0,
        k1 * p - k2 * q * x + k3 * x * x * y - k4 * x,
        k2 * q * x - k3 * x * x * y,
        k2 * q * x,
        k4 * x
    ])


initials = [0.59, 0.56, 0.83, 1.25, 0, 0]  # list of starting values of the variables; first part of the parameters.
coefficients = [0.19, 1.07, 0.85, 0.22]  # list of coefficients for reaction speeds; second part of the parameters.
interval = (0, 200)  # cutoff point in time to stop the simulation at, or None for the default value of 50.
granularity = None  # number of points in time to actually log the values at (not counting t=0),
# or None to let the solver itself decide for us.

if granularity is None:
    evaluations_list = None
else:
    evaluations_list = [np.linspace(interval[0], interval[1], granularity + 1)]

result = simulate(f, [initials], [coefficients], [interval], evaluations_list)
f, initials_list, coefficients_list, interval_list, evaluations_list, solution_list = result

fig, axs = plt.subplots(3, 2, layout='constrained')
for n in range(len(solution_list)):
    initials_out = initials_list[n]
    coefficients_out = coefficients_list[n]
    solution_out = solution_list[n]
    t = solution_out.t
    p, q, x, y, r, s = solution_out.y
    p0, q0, x0, y0, r0, s0 = initials_out
    k1, k2, k3, k4 = coefficients_out
    the_label = f"u0=({p0}, {q0}, {x0}, {y0}, {r0}, {s0});\n{k1=}; {k2=}; {k3=}; {k4=}"
    axs[0, 0].plot(t, p, label=the_label)
    axs[0, 1].plot(t, q, label=the_label)
    axs[1, 0].plot(t, x, label=the_label)
    axs[1, 1].plot(t, y, label=the_label)
    axs[2, 0].plot(t, r, label=the_label)
    axs[2, 1].plot(t, s, label=the_label)
for n in range(len(solution_list)):
    initials_out = initials_list[n]
    coefficients_out = coefficients_list[n]
    solution_out = solution_list[n]
    t = solution_out.t
    p, q, x, y, r, s = solution_out.y
    p0, q0, x0, y0, r0, s0 = initials_out
    k1, k2, k3, k4 = coefficients_out
    the_label = f"u0=({p0}, {q0}, {x0}, {y0}, {r0}, {s0});\n{k1=}; {k2=}; {k3=}; {k4=}"
    axs[1, 0].plot(t, y, label=the_label, alpha=0.2)  # for comparison
    axs[1, 1].plot(t, x, label=the_label, alpha=0.2)  # for comparison
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
axs[1, 0].set_ylim(min(x_min, y_min), max(x_max, y_max))
axs[1, 1].set_ylim(min(x_min, y_min), max(x_max, y_max))
axs[2, 0].set_ylim(r_min, r_max)
axs[2, 1].set_ylim(s_min, s_max)
print("graph success")
handles, labels = axs[2, 0].get_legend_handles_labels()
fig.legend(handles, labels, mode='expand', loc='outside lower center')
fig.suptitle("A graph of the unreduced base model with constant P and Q.")
plt.show()
