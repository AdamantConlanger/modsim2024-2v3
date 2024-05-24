import scipy.integrate as scint
import numpy as np

from typing import Optional, Iterable
flint = float | int
interval_t = tuple[flint, flint] | Iterable[flint]


def execute_simulation(f,
                       u0: Iterable[flint],
                       coeffs: Iterable[flint],
                       time_bound: Optional[flint] = 50,
                       offsets: Optional[str | flint | Iterable[Optional[str | flint]]] = 1,
                       offsets_intervals: Optional[interval_t | Iterable[interval_t]] = (0.75, 1.25),
                       offsets_scale: Optional[int] = 100,
                       time_resolution: Optional[int | str] = "auto",
                       iter_total: Optional[int] = 1):
    rng = np.random.default_rng()

    # declare some lists
    u0_list = []
    coeffs_list = []
    sol_list = []

    # declare the total number of parameters
    total_length = len(u0) + len(coeffs)

    # handle default values and declare total_length.
    # also make sure this function is effectively pass-by-value.
    time_bound_default = 50
    offsets_default = 1
    individual_offset_default = 1
    offsets_intervals_default = (0.75, 1.25)
    offsets_scale_default = 100
    time_resolution_default = "auto"
    iter_total_default = 1

    u0_base = u0
    coeffs_base = coeffs
    time_bound_base = time_bound_default if time_bound is None else time_bound
    offsets_base = offsets_default if offsets is None else offsets
    offsets_intervals_base = offsets_intervals_default if offsets_intervals is None else offsets_intervals
    offsets_scale_base = offsets_scale_default if offsets_scale is None else offsets_scale
    time_resolution_base = time_resolution_default if time_resolution is None else time_resolution
    iter_total_base = iter_total_default if iter_total is None else iter_total

    u0_base = [item for item in u0_base].copy()
    coeffs_base = [item for item in coeffs_base].copy()
    offsets_base = ["random" for i in range(total_length)] if offsets_base == "random" else offsets_base
    if isinstance(offsets_base, str):
        raise ValueError("`offsets` must be either \"random\" or a non-string, which is violated here.")
    offsets_base = [offsets_base for i in range(total_length)] if isinstance(offsets_base, flint) else offsets_base
    offsets_base = [individual_offset_default if item is None else item for item in offsets_base]
    if any([isinstance(item, str) and item != "random" for item in offsets_base]):
        raise ValueError("entries of `offsets` must be either \"random\" or non-strings, which is violated here.")
    if len(offsets_intervals_base) == 2 and isinstance(offsets_intervals_base[0], flint):
        offsets_intervals_base = [[item for item in offsets_intervals_base] for i in range(total_length)]
    offsets_intervals_base = [tuple([value for value in item]) for item in offsets_intervals_base]
    if any([(not isinstance(item, tuple)) or len(item) != 2 for item in offsets_intervals_base]):
        raise ValueError("either `offsets_intervals` or an entry thereof is not a 2-tuple, which is not allowed.")
    if any([any([not isinstance(value, flint) for value in item]) for item in offsets_intervals_base]):
        raise ValueError("`offsets_intervals` does not describe one or more numeric intervals, which is not allowed.")
    if isinstance(time_resolution_base, str) and time_resolution_base != "auto":
        raise ValueError("`time_resolution` must be either \"auto\" or a non-string, which is violated here.")

    # run simulations
    for n in range(iter_total_base):
        # randomize starting values and parameters just a little bit if the user wants that
        transposable = np.array(offsets_intervals_base) * offsets_scale_base
        offsets_intervals_low, offsets_intervals_high = np.transpose(transposable)
        offsets_random = rng.integers(offsets_intervals_low, offsets_intervals_high, total_length) / offsets_scale
        modified_offsets = [offsets_random[i] if offsets_base[i] == "random" else offsets_base[i]
                            for i in range(total_length)]

        modified_u0 = list(np.array(u0_base) * np.array(modified_offsets[0:len(u0_base)]))
        print(f"{modified_u0=}")
        modified_coeffs = list(np.array(coeffs_base) * np.array(modified_offsets[len(u0_base):total_length]))
        print(f"{modified_coeffs=}")
        u0_list.append(modified_u0.copy())
        coeffs_list.append(modified_coeffs.copy())

        trb = time_resolution_base

        t_span_input = (0, time_bound_base)
        y0_input = modified_u0
        args_input = (modified_coeffs,)
        method_input = 'Radau'
        t_eval_input = None if trb == "auto" else np.linspace(0, time_bound_base, trb + 1)

        # try to solve the ODE problem
        try:
            sol = scint.solve_ivp(fun=f,
                                  t_span=t_span_input,
                                  y0=y0_input,
                                  method=method_input,
                                  t_eval=t_eval_input,
                                  args=args_input)
            sol_list.append(sol)
        except Exception:
            print("fail")
            return (False, u0_base, coeffs_base, time_bound_base, time_resolution_base, iter_total_base, u0_list, coeffs_list, sol_list)
        print("success")

    # return the results
    return (True, u0_base, coeffs_base, time_bound_base, time_resolution_base, iter_total_base, u0_list, coeffs_list, sol_list)
