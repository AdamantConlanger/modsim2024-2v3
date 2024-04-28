import scipy.integrate as scint
import numpy as np
import matplotlib.pyplot as plt

from typing import Optional, Iterable
flint = float | int
interval_t = tuple[flint, flint] | Iterable[flint]

def get_variables_from_u(u):
    return [[elem[i] for elem in u] for i in range(len(u[0]))]


def execute_simulation(f,
                       u0: Optional[Iterable[flint]] = None,
                       coeffs: Optional[Iterable[flint]] = None,
                       time_bound: Optional[flint] = 25,
                       offsets: Optional[str | flint | Iterable[Optional[str | flint]]] = "random",
                       offsets_intervals: Optional[interval_t | Iterable[interval_t]] = (0.75, 1.25),
                       offsets_scale: Optional[int] = 1000,
                       iter_total: Optional[int] = 1):
    rng = np.random.default_rng()

    # declare some lists
    u0_list = []
    coeffs_list = []
    sol_list = []

    # handle default values and declare total_length.
    # also make sure this function is effectively pass-by-value.
    u0_default = [1., 1., 0.5, 0.5, 1., 1.]
    coeffs_default = [0.2, 1., 1., 0.2]

    u0_base = u0_default if u0 is None else u0
    coeffs_base = coeffs_default if coeffs is None else coeffs

    total_length = len(u0) + len(coeffs)

    time_bound_default = 25
    offsets_default = "random"
    offsets_intervals_default = (0.75, 1.25)
    offsets_scale_default = 1000
    iter_total_default = 1

    time_bound_base = time_bound_default if time_bound is None else time_bound
    offsets_base = offsets_default if offsets is None else offsets
    offsets_intervals_base = offsets_intervals_default if offsets_intervals is None else offsets_intervals
    offsets_scale_base = offsets_scale_default if offsets_scale is None else offsets_scale
    iter_total_base = iter_total_default if iter_total is None else iter_total

    u0_base = [item for item in u0_base].copy()
    coeffs_base = [item for item in coeffs_base].copy()
    offsets_base = ["random" for i in range(total_length)] if offsets_base == "random" else offsets_base
    if isinstance(offsets_base, str):
        raise ValueError("`offsets` must be either \"random\" or a non-string, which is violated here.")
    offsets_base = [offsets_base for i in range(total_length)] if isinstance(offsets_base, flint) else offsets_base
    offsets_base = [None if item == "random" else item for item in offsets_base]
    if any([isinstance(item, str) for item in offsets_base]):
        raise ValueError("entries of `offsets` must be either \"random\" or non-strings, which is violated here.")
    if len(offsets_intervals_base) == 2 and isinstance(offsets_intervals_base[0], flint):
        offsets_intervals_base = [[item for item in offsets_intervals_base] for i in range(total_length)]
    offsets_intervals_base = [tuple([value for value in item]) for item in offsets_intervals_base]
    if any([(not isinstance(item, tuple)) or len(item) != 2 for item in offsets_intervals_base]):
        raise ValueError("either `offsets_intervals` or an entry thereof is not a 2-tuple, which is not allowed.")
    if any([any([not isinstance(value, flint) for value in item]) for item in offsets_intervals_base]):
        raise ValueError("`offsets_intervals` does not describe one or more numeric intervals, which is not allowed.")

    # run simulations
    for n in range(iter_total_base):
        # randomize starting values and parameters just a little bit if the user wants that
        transposable = np.array(offsets_intervals_base) * offsets_scale_base
        offsets_intervals_low, offsets_intervals_high = np.transpose(transposable)
        offsets_random = rng.integers(offsets_intervals_low, offsets_intervals_high, total_length) / offsets_scale
        modified_offsets = [offsets_random[i] if offsets_base[i] is None else offsets_base[i]
                            for i in range(total_length)]

        modified_u0 = list(np.array(u0_base) * np.array(modified_offsets[0:len(u0_base)]))
        print(f"{modified_u0=}")
        modified_coeffs = list(np.array(coeffs_base) * np.array(modified_offsets[len(u0_base):total_length]))
        print(f"{modified_coeffs=}")
        u0_list.append(modified_u0.copy())
        coeffs_list.append(modified_coeffs.copy())

        t_span_input = (0, time_bound_base)
        y0_input = modified_u0
        args_input = (modified_coeffs,)
        method_input = 'Radau'

        # try to solve the ODE problem
        try:
            sol = scint.solve_ivp(fun=f,
                                  t_span=t_span_input,
                                  y0=y0_input,
                                  args=args_input,
                                  method=method_input)
            sol_list.append(sol)
        except:
            print("fail")
            return (False, u0_base, coeffs_base, time_bound_base, iter_total_base, u0_list, coeffs_list, sol_list)
        print("success")

    # return the results
    return (True, u0_base, coeffs_base, time_bound_base, iter_total_base, u0_list, coeffs_list, sol_list)


def visualize_simulation(u0_base, coeffs_base, time_bound_base, iter_total_base, u0_list, coeffs_list, sol_list):
    fig, axs = plt.subplots(3, 2, layout='constrained')
    try:
        for n in range(len(sol_list)):
            u0 = u0_list[n]
            coeffs = coeffs_list[n]
            sol = sol_list[n]
            t = sol.t
            p, q, x, y, r, s = sol.y
            p0, q0, x0, y0, r0, s0 = u0
            k1, k2, k3, k4 = coeffs
            the_label = f"u0=({p0}, {q0}, {x0}, {y0}, {r0}, {s0});\n{k1=}; {k2=}; {k3=}; {k4=}"
            axs[0, 0].plot(t, p, label=the_label)
            axs[0, 1].plot(t, q, label=the_label)
            axs[1, 0].plot(t, x, label=the_label)
            axs[1, 1].plot(t, y, label=the_label)
            axs[2, 0].plot(t, r, label=the_label)
            axs[2, 1].plot(t, s, label=the_label)
        for n in range(len(sol_list)):
            u0 = u0_list[n]
            coeffs = coeffs_list[n]
            sol = sol_list[n]
            t = sol.t
            p, q, x, y, r, s = sol.y
            p0, q0, x0, y0, r0, s0 = u0
            k1, k2, k3, k4 = coeffs
            the_label = f"u0=({p0}, {q0}, {x0}, {y0}, {r0}, {s0});\n{k1=}; {k2=}; {k3=}; {k4=}"
            axs[1, 0].plot(t, y, label=the_label, alpha=0.2)  # for comparison
            axs[1, 1].plot(t, x, label=the_label, alpha=0.2)  # for comparison
            # TODO: make the color of these lines appear in the legend too
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
        axs[0, 0].set_ylim(0, 2 * u0_base[0])
        axs[0, 1].set_ylim(0, 2 * u0_base[1])
        axs[1, 0].set_ylim(min(x_min, y_min), max(x_max, y_max))
        axs[1, 1].set_ylim(min(x_min, y_min), max(x_max, y_max))
        axs[2, 0].set_ylim(r_min, r_max)
        axs[2, 1].set_ylim(s_min, s_max)
    except:
        print("graph fail")
    print("graph success")
    # for ax in axs.flat:
    #    ax.label_outer()
    handles, labels = axs[2, 0].get_legend_handles_labels()
    fig.legend(handles, labels, mode='expand', loc='outside lower center')
    fig.suptitle("A graph of the base model with constant P and Q.")
    plt.show()


def f(t, u, coeffs):
    p, q, x, y, r, s = u
    k1, k2, k3, k4 = coeffs
    return np.array([
        0,
        0,
        k1 * p - k2 * q * x + k3 * x**2 * y - k4 * x,
        k2 * q * x - k3 * x**2 * y,
        k2 * q * x,
        k4 * x
    ])


# u0 = [0.5, 0.5, 1, 1, 1, 1]
# coeffs = [0.2, 1, 1, 0.2]
u0 = [0.59, 0.56, 0.83, 1.25, 0, 0]
coeffs = [0.19, 1.07, 0.85, 0.22]
time_bound = 200
offsets = [1 for i in range(0, 10)]
# offsets = "random"
iter_total = None
# iter_total = 4

result = execute_simulation(f, u0=u0, coeffs=coeffs, time_bound=time_bound, offsets=offsets, iter_total=iter_total)
success, u0_base, coeffs_base, time_bound_base, iter_total_base, u0_list, coeffs_list, sol_list = result

if success:
    visualize_simulation(u0_base, coeffs_base, time_bound_base, iter_total_base, u0_list, coeffs_list, sol_list)
else:
    print("simulation failed. Stopping and printing failed settings.")
    print(u0_list, "\n", coeffs_list)
