import math
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt


# takes in the output of a TSQ simulation and determines whether or not it collapsed,
# and what the amplitude (on both sides) and period are if not.
# collapse is defined as staying within the tolerance limits of the equilibrium
# for the last 25% of the simulation.
# if, however, the period of time in which x and y are within tolerances at the end
# is not the longest span of time that x and y have stayed within tolerances for in the entire sim,
# we deem it uncertain whether or not they have actually collapsed.
# note: the actual algorithm for determining this may be changed without it being mentioned here.
# tolerance is a tuple of the x-tolerance and y-tolerance
# if relative tolerances are desired, equi * tolerance is used as true tolerance limit radius.
# Else, tolerance is used.


def determine_and_plot_metrics(coefficients_list, solution_list, separate_modified_variables, coeffs_extents, tolerance, use_relative=False, plot_periods=True, plot_upper_x_amps=False, plot_lower_x_amps=False, plot_upper_y_amps=False, plot_lower_y_amps=False, plot_collapse_moments=False):
    metrics_list = []
    for n in range(len(solution_list)):
        current_metrics = dict()
        coefficients_out = coefficients_list[n]
        solution_out = solution_list[n]
        t = solution_out.t
        x, y = solution_out.y
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

            # free some memory:
            solution_list[n] = None

        metrics_list.append(current_metrics)




        # freeing memory
        current_metrics = None
        coefficients_out = None
        solution_out = None
        t = None
        x = None
        y = None
        mu1 = None
        x_equi = None
        y_equi = None
        tol_x = None
        tol_y = None
        tolerability_x = None
        tolerability_y = None
        tolerability = None
        late_index = None
        seems_collapsed = None
        tolerability_bounds_reverse_indices = None
        tolerability_span_lengths = None
        lower_index = None
        upper_index = None
        seems_collapsed = None
        is_collapsed = None
        is_collapse_certain = None
        collapse_moment = None
        use_exits_x = None
        use_exits_y = None
        x_events = None
        y_events = None
        use_falling_x = None
        use_falling_y = None
        x_usable_events = None
        y_usable_events = None
        use_events_x = None
        use_events_y = None
        x_period_estimates = None
        x_weighted_period_estimates = None
        x_weights = None
        x_period = None
        y_period_estimates = None
        y_weighted_period_estimates = None
        y_weights = None
        y_period = None
        estimated_period = None
        x_maxima_indices = None
        x_minima_indices = None
        y_maxima_indices = None
        y_minima_indices = None
        is_x_falling = None
        is_y_falling = None
        is_x_rising = None
        is_y_rising = None
        x_has_maximum = None
        x_has_minimum = None
        y_has_maximum = None
        y_has_minimum = None
        x_maximum = None
        x_minimum = None
        y_maximum = None
        y_minimum = None
        x_upper_amplitude = None
        x_lower_amplitude = None
        y_upper_amplitude = None
        y_lower_amplitude = None

    metrics = metrics_list


    # freeing memory
    solution_list = None

    # coefficients_to_indices = dict()
    mu1_list = separate_modified_variables[2]
    mu2_list = separate_modified_variables[3]
    mu1_count = len(mu1_list)
    mu2_count = len(mu2_list)

    if plot_periods:
        periods = [item["estimated_period"] if "estimated_period" in item else 0 if item["is_collapsed"] else np.nan for item in metrics]
    elif plot_upper_x_amps:
        upper_x_amps = [item["x_upper_amplitude"]
                        if "x_upper_amplitude" in item else 0 if item["is_collapsed"] else np.nan for item in metrics]
    elif plot_lower_x_amps:
        lower_x_amps = [item["x_lower_amplitude"]
                        if "x_lower_amplitude" in item else 0 if item["is_collapsed"] else np.nan for item in metrics]
    elif plot_upper_y_amps:
        upper_y_amps = [item["y_upper_amplitude"]
                        if "y_upper_amplitude" in item else 0 if item["is_collapsed"] else np.nan for item in metrics]
    elif plot_lower_y_amps:
        lower_y_amps = [item["y_lower_amplitude"]
                        if "y_lower_amplitude" in item else 0 if item["is_collapsed"] else np.nan for item in metrics]
    elif plot_collapse_moments:
        collapse_moments = [item["collapse_moment"] if "collapse_moment" in item else np.nan for item in metrics]
    else:
        return


    # freeing memory
    metrics = None

    mu1_min, mu1_max, mu2_min, mu2_max = coeffs_extents
    desired_range_x, desired_range_y = np.meshgrid(np.linspace(mu1_min, mu1_max, mu1_count),
                                                np.linspace(mu2_min, mu2_max, mu2_count), indexing='ij')
    print(f"{mu1_list=}")
    print(f"{mu2_list=}")
    print(f"{desired_range_x=}")
    print(f"{desired_range_y=}")
    print(f"{coefficients_list=}")
    desired_range = (desired_range_x, desired_range_y)

    # https://stackoverflow.com/a/14140554/18375328,
    # which is superseded by https://docs.scipy.org/doc/scipy/tutorial/interpolate/ND_unstructured.html

    # TODO: make this dependent on which stuff needs to be plotted

    plottable_data = griddata(coefficients_list, periods, desired_range, method="nearest")


    fig, axs = plt.subplots(1, 1, layout='constrained')

    use_smooth_block_interpolation = False
    use_smooth_lossy_interpolation = False

    if use_smooth_block_interpolation:
        im = axs.imshow(plottable_data.T, extent=coeffs_extents, interpolation="spline16", origin='lower', cmap="viridis")
    elif use_smooth_lossy_interpolation:
        im = axs.imshow(plottable_data.T, extent=coeffs_extents, interpolation="bicubic", origin='lower', cmap="viridis")
    else:
        im = axs.imshow(plottable_data.T, extent=coeffs_extents, interpolation="none", origin='lower', cmap="viridis")
    axs.set_ylabel("mu2")
    axs.set_xlabel("mu1")
    fig.suptitle("A graph of the period of oscillation for the truncated simple (reduced) quoduct model.\nx0=1, y0=1.")
    fig.colorbar(im, ax=axs, label='period of oscillation (0 if n/a)')
    plt.show()

            
                    






            
            


            
            


