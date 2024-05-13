import math
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import scipy.integrate as scint


def simulate_and_plot_metrics(f, initials_list, coefficients_list, interval_list, tolerance=(10**-4, 10**-4), evaluations_list=None, use_relative=False, plot_periods=True, plot_upper_x_amps=False, plot_lower_x_amps=False, plot_upper_y_amps=False, plot_lower_y_amps=False, plot_collapse_moments=False):
    mu1_max = max([item[0] for item in coefficients_list])
    mu1_min = min([item[0] for item in coefficients_list])
    mu2_max = max([item[1] for item in coefficients_list])
    mu2_min = min([item[1] for item in coefficients_list])
    mu1_len = len(set([item[0] for item in coefficients_list]))
    mu2_len = len(set([item[1] for item in coefficients_list]))
    extent = (mu1_min, mu1_max, mu2_min, mu2_max)
    

    metrics_list = []

    for iteration in range(len(initials_list)):
        t_eval = None if evaluations_list is None else evaluations_list[iteration]

        sol = scint.solve_ivp(fun=f,
                              t_span=interval_list[iteration],
                              y0=initials_list[iteration],
                              method="Radau",
                              t_eval=t_eval,
                              args=(coefficients_list[iteration],))

        current_metrics = dict()
        coefficients_out = coefficients_list[iteration]
        t = sol.t
        x, y = sol.y
        mu1 = coefficients_out[0]
        x_equi = mu1
        y_equi = 1/mu1
        tol_x = tolerance[0] * x_equi if use_relative else tolerance[0]
        tol_y = tolerance[1] * y_equi if use_relative else tolerance[1]

        # state that everything's unknown if there are fewer than 100 evaluations.
        if len(t) < 100:
            current_metrics["unknown"] = True
            metrics_list.append(current_metrics)
            continue
        else:
            current_metrics["unknown"] = False

        # determine when x and y are "near enough" to the equilibria
        tolerability_x = [abs(item - x_equi) <= tol_x for item in x]
        tolerability_y = [abs(item - y_equi) <= tol_y for item in y]
        tolerability = [tolerability_x[index] and tolerability_y[index] for index in range(len(t))]

        # determine whether collapse: first check whether the last 25% is within tolerance
        late_index = math.floor(len(t) * 3 / 4)
        seems_collapsed = any(tolerability[late_index:])

        # if it collapses, check whether that's the longest span and deem it certain or uncertain
        if seems_collapsed:
            # determine the boundaries of the regions where x and y are tolerable
            tolerability_bounds_reverse_indices = [0]
            for index in range(1, len(t)):
                if tolerability[-index] != tolerability[-(index+1)]:
                    tolerability_bounds_reverse_indices.append(index)
            if len(tolerability_bounds_reverse_indices) % 2 != 0:
                # if x[0] and y[0] are within tolerances, then that's a boundary point.
                tolerability_bounds_reverse_indices.append(len(t))
            # determine the lengths of these regions
            tolerability_span_lengths = []
            for index in range(len(tolerability_bounds_reverse_indices)):
                if index % 2 == 0:
                    lower_index = tolerability_bounds_reverse_indices[index]
                    upper_index = tolerability_bounds_reverse_indices[index + 1]
                    tolerability_span_lengths.append(upper_index - lower_index)
            # if the final span is strictly the longest, we're pretty sure it collapses. Otherwise, we aren't.
            if len(tolerability_span_lengths) == 1:
                seems_collapsed = True
                is_collapsed = True
                is_collapse_certain = True
            elif tolerability_span_lengths[0] > 1.5 * max(tolerability_span_lengths[1:]):
                seems_collapsed = True
                is_collapsed = True
                is_collapse_certain = True
            elif tolerability_span_lengths[0] > max(tolerability_span_lengths[1:]):
                seems_collapsed = True
                is_collapsed = True
                is_collapse_certain = False
            else:
                seems_collapsed = True
                is_collapsed = False
                is_collapse_certain = False
            collapse_moment = len(t) - tolerability_bounds_reverse_indices[1]

            # add metrics to current metrics
            current_metrics["collapse_moment"] = collapse_moment

        else:
            # if it doesn't seem collapsed, it also isn't collapsed
            is_collapsed = False
            is_collapse_certain = False

        # add metrics to current metrics
        current_metrics["seems_collapsed"] = seems_collapsed
        current_metrics["is_collapsed"] = is_collapsed
        current_metrics["is_collapse_certain"] = is_collapse_certain

        # seems_collapsed: does it stay within tolerances for the last part of the sim?
        # is_collapsed: is this the first time it's stayed within tolerances for as long as it has?
        # is_collapse_certain: is it staying within tolerances for considerably longer than ever before?
        # collapse_moment: first point within tolerances after which it stays so for the rest of the sim.

        # if it doesn't collapse, we determine periods and amplitudes
        if not is_collapsed:
            # we need to determine all local maxima and minima and "characterisic" points for calculating periods,
            # both for x and for y, separately.

            # if x is tolerable at the very very end of the sim, we'll check for points where x becomes tolerable.
            # otherwise, we look for points where x leaves the range of tolerability. Same for y.
            use_exits_x = not tolerability_x[-1]
            use_exits_y = not tolerability_y[-1]

            # first, we determine all the places x and y leave/enter their ranges of tolerability (separately),
            # or full-on cross the equilibrium.
            x_events = []
            y_events = []
            for index in range(len(t) - 1):
                if x[index] > x_equi and x[index + 1] < x_equi:
                    x_events.append(index)
                    continue
                if x[index] < x_equi and x[index + 1] > x_equi:
                    x_events.append(index)
                    continue
                if use_exits_x and tolerability_x[index] and not tolerability_x[index + 1]:
                    x_events.append(index)
                if not use_exits_x and not tolerability_x[index] and tolerability_x[index + 1]:
                    x_events.append(index)

            for index in range(len(t) - 1):
                if y[index] > y_equi and y[index + 1] < y_equi:
                    y_events.append(index)
                    continue
                if y[index] < y_equi and y[index + 1] > y_equi:
                    y_events.append(index)
                    continue
                if use_exits_y and tolerability_y[index] and not tolerability_y[index + 1]:
                    y_events.append(index)
                if not use_exits_y and not tolerability_y[index] and tolerability_y[index + 1]:
                    y_events.append(index)

            # first, we just give up if either x or y never reaches or passes through the equilibrium
            if len(x_events) == 0 or len(y_events) == 0:
                current_metrics["has_usable_oscillations"] = False
                metrics_list.append(current_metrics)
                continue

            # next, we determine the monotonicity of x and y at the respective last of these points.
            use_falling_x = x[x_events[-1] + 1] > x[x_events[-1]]
            use_falling_y = y[y_events[-1] + 1] > y[y_events[-1]]

            # we only consider the events that have the correct monotonicities
            x_usable_events = []
            y_usable_events = []
            for index in x_events:
                if use_falling_x and x[index + 1] > x[index]:
                    x_usable_events.append(index)
                if not use_falling_x and x[index + 1] < x[index]:
                    x_usable_events.append(index)
            for index in y_events:
                if use_falling_y and y[index + 1] > y[index]:
                    y_usable_events.append(index)
                if not use_falling_y and y[index + 1] < y[index]:
                    y_usable_events.append(index)

            # we now check whether there are at least two such events
            use_events_x = len(x_usable_events) > 1
            use_events_y = len(y_usable_events) > 1

            # we now estimate the period based on a weighted average of the time between these events,
            # with later events being considered better estimates than earlier ones.
            if use_events_x:
                x_period_estimates = []
                for index in range(1, len(x_usable_events)):
                    x_period_estimates.append(t[x_usable_events[index]] - t[x_usable_events[index - 1]])
                x_weighted_period_estimates = []
                x_weights = []
                for index in range(len(x_period_estimates)):
                    x_weighted_period_estimates.append(x_period_estimates[index] / (len(x_period_estimates) - index))
                    x_weights.append(1 / (len(x_period_estimates) - index))
                x_period = sum(x_weighted_period_estimates) / sum(x_weights)
            if use_events_y:
                y_period_estimates = []
                for index in range(1, len(y_usable_events)):
                    y_period_estimates.append(t[y_usable_events[index]] - t[y_usable_events[index - 1]])
                y_weighted_period_estimates = []
                y_weights = []
                for index in range(len(y_period_estimates)):
                    y_weighted_period_estimates.append(y_period_estimates[index] / (len(y_period_estimates) - index))
                    y_weights.append(1 / (len(y_period_estimates) - index))
                y_period = sum(y_weighted_period_estimates) / sum(y_weights)

            # we now take the mean of of the x estimate and the y estimate for the period
            if use_events_x and use_events_y:
                estimated_period = (x_period + y_period) / 2
            elif use_events_x:
                estimated_period = x_period
            elif use_events_y:
                estimated_period = y_period
            else:
                # if neither x nor y is usable, just give up
                current_metrics["has_usable_oscillations"] = False
                metrics_list.append(current_metrics)
                continue
            current_metrics["has_usable_oscillations"] = True

            # add the period to the metrics
            current_metrics["estimated_period"] = estimated_period

            # next, we calculate the local minima and maxima of x and y.
            x_maxima_indices = []
            x_minima_indices = []
            y_maxima_indices = []
            y_minima_indices = []
            is_x_falling = x[1] <= x[0]
            is_y_falling = y[1] <= y[0]
            is_x_rising = x[1] >= x[0]
            is_y_rising = y[1] >= y[0]
            for index in range(1, len(t)):
                if x[index] > x[index - 1]:
                    if is_x_falling:
                        x_minima_indices.append(index - 1)
                        is_x_falling = False
                    is_x_rising = True
                if x[index] < x[index - 1]:
                    if is_x_rising:
                        x_maxima_indices.append(index - 1)
                        is_x_rising = False
                    is_x_falling = True
                if y[index] > y[index - 1]:
                    if is_y_falling:
                        y_minima_indices.append(index - 1)
                        is_y_falling = False
                    is_y_rising = True
                if y[index] < y[index - 1]:
                    if is_y_rising:
                        y_maxima_indices.append(index - 1)
                        is_y_rising = False
                    is_y_falling = True

            # if there are no minima and/or maxima, give up on that department
            x_has_maximum = len(x_maxima_indices) > 0
            x_has_minimum = len(x_minima_indices) > 0
            y_has_maximum = len(y_maxima_indices) > 0
            y_has_minimum = len(y_minima_indices) > 0

            current_metrics["x_has_maximum"] = x_has_maximum
            current_metrics["x_has_minimum"] = x_has_minimum
            current_metrics["y_has_maximum"] = y_has_maximum
            current_metrics["y_has_minimum"] = y_has_minimum

            # if possible and reasonable, we disregard the first half of the local extrema
            if x_has_maximum:
                x_maxima_indices = x_maxima_indices[math.floor(len(x_maxima_indices) / 3):]
            if x_has_minimum:
                x_minima_indices = x_minima_indices[math.floor(len(x_minima_indices) / 3):]
            if y_has_maximum:
                y_maxima_indices = y_maxima_indices[math.floor(len(y_maxima_indices) / 3):]
            if y_has_minimum:
                y_minima_indices = y_minima_indices[math.floor(len(y_minima_indices) / 3):]

            # we now take the largest maxima and smallest minima as our true maxima
            if x_has_maximum:
                x_maximum = max([x[index] for index in x_maxima_indices])
            if x_has_minimum:
                x_minimum = max([x[index] for index in x_minima_indices])
            if y_has_maximum:
                y_maximum = max([y[index] for index in y_maxima_indices])
            if y_has_minimum:
                y_minimum = max([y[index] for index in y_minima_indices])

            # finally, we calculate the amplitudes
            if x_has_maximum:
                x_upper_amplitude = x_maximum - x_equi
            if x_has_minimum:
                x_lower_amplitude = x_equi - x_minimum
            if y_has_maximum:
                y_upper_amplitude = y_maximum - y_equi
            if y_has_minimum:
                y_lower_amplitude = y_equi - y_minimum

            # and we add the amplitudes to the metrics
            current_metrics["x_upper_amplitude"] = x_upper_amplitude if x_has_maximum else None
            current_metrics["x_lower_amplitude"] = x_lower_amplitude if x_has_minimum else None
            current_metrics["y_upper_amplitude"] = y_upper_amplitude if y_has_maximum else None
            current_metrics["y_lower_amplitude"] = y_lower_amplitude if y_has_minimum else None

        metrics_list.append(current_metrics)

    desired_range_x, desired_range_y = np.meshgrid(np.linspace(mu1_min, mu1_max, mu1_len),
                                                   np.linspace(mu2_min, mu2_max, mu2_len), indexing='ij')
    desired_range = (desired_range_x, desired_range_y)
    print(desired_range)
    print(coefficients_list)

    metrics_range = range(len(metrics_list))
    are_collapsed = [item["is_collapsed"] if "is_collapsed" in item else False for item in metrics_list]


    

    if plot_periods:
        have_period = ["estimated_period" in item for item in metrics_list]
        periods = [metrics_list[index]["estimated_period"] if have_period[index] else np.nan for index in metrics_range]
        periods = [0 if are_collapsed[index] else periods[index] for index in metrics_range]
        periods_data = griddata(coefficients_list, periods, desired_range, method="nearest")

        fig, axs = plt.subplots(1, 1, layout='constrained')
        dx, = np.diff(extent[:2])/(periods_data.T.shape[1]-1)
        dy, = -np.diff(extent[2:])/(periods_data.T.shape[0]-1)
        shifted_extent = [extent[0]-dx/2, extent[1]+dx/2, extent[2]+dy/2, extent[3]-dy/2]
        im = axs.imshow(periods_data.T, interpolation="none", origin='lower',
                        cmap="viridis", extent=shifted_extent, aspect='auto')
        axs.set_xticks(np.arange(extent[0], extent[1]+dx, dx))
        axs.set_yticks(np.arange(extent[3], extent[2]+dy, dy))
        axs.set_ylabel("mu2")
        axs.set_xlabel("mu1")
        fig.suptitle("A graph of the period of oscillation\nfor the truncated simple (reduced) quoduct model,\nwith x0=1, y0=1.")
        fig.colorbar(im, ax=axs, label='period of oscillation (0 if converges)')
        plt.show()
    if plot_upper_x_amps:
        have_x_upper_amps = ["x_upper_amplitude" in item for item in metrics_list]
        x_upper_amps = [metrics_list[index]["x_upper_amplitude"]
                   if have_x_upper_amps[index] else np.nan for index in metrics_range]
        x_upper_amps = [0 if are_collapsed[index] else x_upper_amps[index] for index in metrics_range]
        x_upper_amps_data = griddata(coefficients_list, x_upper_amps, desired_range, method="nearest")

        fig, axs = plt.subplots(1, 1, layout='constrained')
        dx, = np.diff(extent[:2])/(x_upper_amps_data.T.shape[1]-1)
        dy, = -np.diff(extent[2:])/(x_upper_amps_data.T.shape[0]-1)
        shifted_extent = [extent[0]-dx/2, extent[1]+dx/2, extent[2]+dy/2, extent[3]-dy/2]
        im = axs.imshow(x_upper_amps_data.T, interpolation="none", origin='lower',
                        cmap="viridis", extent=shifted_extent, aspect='auto')
        axs.set_xticks(np.arange(extent[0], extent[1]+dx, dx))
        axs.set_yticks(np.arange(extent[3], extent[2]+dy, dy))
        axs.set_ylabel("mu2")
        axs.set_xlabel("mu1")
        fig.suptitle("A graph of the upper x wave amplitude\nfor the truncated simple (reduced) quoduct model,\nwith x0=1, y0=1.")
        fig.colorbar(im, ax=axs, label='upper x amplitude (0 if converges)')
        plt.show()
    if plot_lower_x_amps:
        have_x_lower_amps = ["x_lower_amplitude" in item for item in metrics_list]
        x_lower_amps = [metrics_list[index]["x_lower_amplitude"]
                        if have_x_lower_amps[index] else np.nan for index in metrics_range]
        x_lower_amps = [0 if are_collapsed[index] else x_lower_amps[index] for index in metrics_range]
        x_lower_amps_data = griddata(coefficients_list, x_lower_amps, desired_range, method="nearest")

        fig, axs = plt.subplots(1, 1, layout='constrained')
        dx, = np.diff(extent[:2])/(x_lower_amps_data.T.shape[1]-1)
        dy, = -np.diff(extent[2:])/(x_lower_amps_data.T.shape[0]-1)
        shifted_extent = [extent[0]-dx/2, extent[1]+dx/2, extent[2]+dy/2, extent[3]-dy/2]
        im = axs.imshow(x_lower_amps_data.T, interpolation="none", origin='lower',
                        cmap="viridis", extent=shifted_extent, aspect='auto')
        axs.set_xticks(np.arange(extent[0], extent[1]+dx, dx))
        axs.set_yticks(np.arange(extent[3], extent[2]+dy, dy))
        axs.set_ylabel("mu2")
        axs.set_xlabel("mu1")
        fig.suptitle("A graph of the lower x wave amplitude\nfor the truncated simple (reduced) quoduct model,\nwith x0=1, y0=1.")
        fig.colorbar(im, ax=axs, label='lower x amplitude (0 if converges)')
        plt.show()
    if plot_upper_y_amps:
        have_y_upper_amps = ["y_upper_amplitude" in item for item in metrics_list]
        y_upper_amps = [metrics_list[index]["y_upper_amplitude"]
                        if have_y_upper_amps[index] else np.nan for index in metrics_range]
        y_upper_amps = [0 if are_collapsed[index] else y_upper_amps[index] for index in metrics_range]
        y_upper_amps_data = griddata(coefficients_list, y_upper_amps, desired_range, method="nearest")

        fig, axs = plt.subplots(1, 1, layout='constrained')
        dx, = np.diff(extent[:2])/(y_upper_amps_data.T.shape[1]-1)
        dy, = -np.diff(extent[2:])/(y_upper_amps_data.T.shape[0]-1)
        shifted_extent = [extent[0]-dx/2, extent[1]+dx/2, extent[2]+dy/2, extent[3]-dy/2]
        im = axs.imshow(y_upper_amps_data.T, interpolation="none", origin='lower',
                        cmap="viridis", extent=shifted_extent, aspect='auto')
        axs.set_xticks(np.arange(extent[0], extent[1]+dx, dx))
        axs.set_yticks(np.arange(extent[3], extent[2]+dy, dy))
        axs.set_ylabel("mu2")
        axs.set_xlabel("mu1")
        fig.suptitle("A graph of the upper y wave amplitude\nfor the truncated simple (reduced) quoduct model,\nwith x0=1, y0=1.")
        fig.colorbar(im, ax=axs, label='upper y amplitude (0 if converges)')
        plt.show()
    if plot_lower_y_amps:
        have_y_lower_amps = ["y_lower_amplitude" in item for item in metrics_list]
        y_lower_amps = [metrics_list[index]["y_lower_amplitude"]
                        if have_y_lower_amps[index] else np.nan for index in metrics_range]
        y_lower_amps = [0 if are_collapsed[index] else y_lower_amps[index] for index in metrics_range]
        y_lower_amps_data = griddata(coefficients_list, y_lower_amps, desired_range, method="nearest")

        fig, axs = plt.subplots(1, 1, layout='constrained')
        dx, = np.diff(extent[:2])/(y_lower_amps_data.T.shape[1]-1)
        dy, = -np.diff(extent[2:])/(y_lower_amps_data.T.shape[0]-1)
        shifted_extent = [extent[0]-dx/2, extent[1]+dx/2, extent[2]+dy/2, extent[3]-dy/2]
        im = axs.imshow(y_lower_amps_data.T, interpolation="none", origin='lower',
                        cmap="viridis", extent=shifted_extent, aspect='auto')
        axs.set_xticks(np.arange(extent[0], extent[1]+dx, dx))
        axs.set_yticks(np.arange(extent[3], extent[2]+dy, dy))
        axs.set_ylabel("mu2")
        axs.set_xlabel("mu1")
        fig.suptitle("A graph of the lower y wave amplitude\nfor the truncated simple (reduced) quoduct model,\nwith x0=1, y0=1.")
        fig.colorbar(im, ax=axs, label='lower y amplitude (0 if converges)')
        plt.show()
    if plot_collapse_moments:
        collapse_moments = [item["collapse_moment"] if "collapse_moment" in item else np.nan for item in metrics_list]
        collapse_moments_data = griddata(coefficients_list, collapse_moments, desired_range, method="nearest")

        fig, axs = plt.subplots(1, 1, layout='constrained')
        dx, = np.diff(extent[:2])/(collapse_moments_data.T.shape[1]-1)
        dy, = -np.diff(extent[2:])/(collapse_moments_data.T.shape[0]-1)
        shifted_extent = [extent[0]-dx/2, extent[1]+dx/2, extent[2]+dy/2, extent[3]-dy/2]
        im = axs.imshow(collapse_moments_data.T, interpolation="none", origin='lower',
                        cmap="viridis", extent=shifted_extent, aspect='auto')
        axs.set_xticks(np.arange(extent[0], extent[1]+dx, dx))
        axs.set_yticks(np.arange(extent[3], extent[2]+dy, dy))
        axs.set_ylabel("mu2")
        axs.set_xlabel("mu1")
        tol_type = "relative" if use_relative else "absolute"
        fig.suptitle(f"A graph of \"the moment\" x and y converge\nfor the truncated simple (reduced) quoduct model,\nwith x0=1, y0=1, {tol_type} tolerances={tolerance}.")
        fig.colorbar(im, ax=axs, label='moment of convergence within tolerances')
        plt.show()
    





