def perform_program(simulate, cartesian_product):
    import scipy.integrate as scint
    from scipy.interpolate import griddata
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    import json

    def simulate_and_plot_metrics(f, initials, coefficients_list, interval, granularity, atol=10**-7, rtol=10**-6, tolerance=10**-4, evaluations=None, use_relative=False, plot_periods=True, plot_x_amps=False, plot_y_amps=False, plot_collapse_moments=False):
        x0, y0 = initials
        alpha_max = max([item[0] for item in coefficients_list])
        alpha_min = min([item[0] for item in coefficients_list])
        beta_max = max([item[1] for item in coefficients_list])
        beta_min = min([item[1] for item in coefficients_list])
        alpha_len = len(set([item[0] for item in coefficients_list]))
        beta_len = len(set([item[1] for item in coefficients_list]))
        extent = (alpha_min, alpha_max, beta_min, beta_max)

        metrics_list = []

        for iteration in range(len(coefficients_list)):
            coefficients_out = coefficients_list[iteration]

            sol = scint.solve_ivp(fun=f,
                                  t_span=interval,
                                  y0=initials,
                                  method="Radau",
                                  t_eval=evaluations,
                                  args=(coefficients_list[iteration],),
                                  rtol=rtol,
                                  atol=atol)
            print(f"finished simulating case {coefficients_out}")

            current_metrics = dict()
            t = sol.t
            x, y = sol.y
            alpha, beta = coefficients_out
            x_equi = alpha / beta
            y_equi = beta / alpha
            tol_x = tolerance * x_equi if use_relative else tolerance
            tol_y = tolerance * y_equi if use_relative else tolerance

            # state that everything's unknown if there are fewer than 100 evaluations.
            if len(t) < 100:
                current_metrics["unknown"] = True
                metrics_list.append(current_metrics)
                print(f"finished case {coefficients_out}")
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
                    print(f"finished case {coefficients_out}")
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
                        x_weighted_period_estimates.append(
                            x_period_estimates[index] / (len(x_period_estimates) - index))
                        x_weights.append(1 / (len(x_period_estimates) - index))
                    x_period = sum(x_weighted_period_estimates) / sum(x_weights)
                if use_events_y:
                    y_period_estimates = []
                    for index in range(1, len(y_usable_events)):
                        y_period_estimates.append(t[y_usable_events[index]] - t[y_usable_events[index - 1]])
                    y_weighted_period_estimates = []
                    y_weights = []
                    for index in range(len(y_period_estimates)):
                        y_weighted_period_estimates.append(
                            y_period_estimates[index] / (len(y_period_estimates) - index))
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
                    print(f"finished case {coefficients_out}")
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
                x_has_amplitude = x_has_maximum and x_has_minimum
                y_has_maximum = len(y_maxima_indices) > 0
                y_has_minimum = len(y_minima_indices) > 0
                y_has_amplitude = y_has_maximum and y_has_minimum

                # if possible and reasonable, we disregard the first half of the local extrema.
                # we then take the largest maxima and smallest minima as our true maxima.
                # and finally, we calculate the amplitudes
                if x_has_amplitude:
                    x_maxima_indices = x_maxima_indices[math.floor(len(x_maxima_indices) / 3):]
                    x_minima_indices = x_minima_indices[math.floor(len(x_minima_indices) / 3):]
                    x_maximum = max([x[index] for index in x_maxima_indices])
                    x_minimum = min([x[index] for index in x_minima_indices])
                    x_amplitude = x_maximum - x_minimum
                if y_has_amplitude:
                    y_maxima_indices = y_maxima_indices[math.floor(len(y_maxima_indices) / 3):]
                    y_minima_indices = y_minima_indices[math.floor(len(y_minima_indices) / 3):]
                    y_maximum = max([y[index] for index in y_maxima_indices])
                    y_minimum = min([y[index] for index in y_minima_indices])
                    y_amplitude = y_maximum - y_minimum

                # and we add the amplitude stuff to the metrics
                current_metrics["x_has_amplitude"] = x_has_amplitude
                current_metrics["y_has_amplitude"] = y_has_amplitude
                current_metrics["x_amplitude"] = x_amplitude if x_has_amplitude else None
                current_metrics["y_amplitude"] = y_amplitude if y_has_amplitude else None

            metrics_list.append(current_metrics)
            print(f"finished case {coefficients_out}")

        with open("str_metrics.json", "w") as the_newfile:
            json.dump(metrics_list, the_newfile)

        desired_range_x, desired_range_y = np.meshgrid(np.linspace(alpha_min, alpha_max, alpha_len),
                                                       np.linspace(beta_min, beta_max, beta_len), indexing='ij')
        desired_range = (desired_range_x, desired_range_y)

        metrics_range = range(len(metrics_list))
        are_collapsed = [item["is_collapsed"] if "is_collapsed" in item else False for item in metrics_list]

        sliver = 10**-6
        # or np.finfo(type(extent[0])).eps

        if plot_periods:
            have_period = ["estimated_period" in item for item in metrics_list]
            periods = [metrics_list[index]["estimated_period"]
                       if have_period[index] else np.nan for index in metrics_range]
            periods = [0 if are_collapsed[index] else periods[index] for index in metrics_range]
            periods_data = griddata(coefficients_list, periods, desired_range, method="nearest")

            fig, axs = plt.subplots(1, 1, layout='constrained')
            if np.diff(extent[:2]) <= sliver:
                shifted_extent = [extent[0] - sliver, extent[1] + sliver]
                the_xticks = np.array([extent[0]])
            else:
                dx, = np.diff(extent[:2])/(periods_data.T.shape[1]-1)
                shifted_extent = [extent[0]-dx/2, extent[1]+dx/2]
                x_step = dx
                while 20 * x_step < np.diff(extent[:2]):
                    x_step = 2 * x_step
                the_xticks = np.arange(extent[0], extent[1]+dx, x_step)
            if np.diff(extent[2:]) <= 10**-6:
                shifted_extent.extend([extent[2] - sliver, extent[3] + sliver])
                the_yticks = np.array([extent[2]])
            else:
                dy, = np.diff(extent[2:])/(periods_data.T.shape[0]-1)
                shifted_extent.extend([extent[2]-dy/2, extent[3]+dy/2])
                y_step = dy
                while 20 * y_step < np.diff(extent[2:]):
                    y_step = 2 * y_step
                the_yticks = np.arange(extent[2], extent[3]+dy, y_step)
            im = axs.imshow(periods_data.T, interpolation="none", origin='lower',
                            cmap="viridis", extent=shifted_extent, aspect='auto')
            axs.set_xticks(the_xticks)
            axs.set_yticks(the_yticks)
            axs.set_ylabel("beta")
            axs.set_xlabel("alpha")
            tol_type = "rel" if use_relative else "abs"
            the_title = "A graph of the period of oscillation\nfor the simple truncated reduced model,\n"
            the_title += f"with {x0=}, {y0=}, t_max={interval[1]}, grains={granularity}," + "\n"
            the_title += f"abs sim tol={atol}, rel sim tol={rtol}, {tol_type} analysis tol={tolerance}."
            fig.suptitle(the_title)
            fig.colorbar(im, ax=axs, label='period of oscillation (0 if converges)')
            plt.show()
        if plot_x_amps:
            have_x_amps = ["x_amplitude" in item for item in metrics_list]
            x_amps = [metrics_list[index]["x_amplitude"]
                      if have_x_amps[index] else np.nan for index in metrics_range]
            x_amps = [0 if are_collapsed[index] else x_amps[index] for index in metrics_range]
            x_amps_data = griddata(coefficients_list, x_amps, desired_range, method="nearest")

            fig, axs = plt.subplots(1, 1, layout='constrained')
            if np.diff(extent[:2]) <= sliver:
                shifted_extent = [extent[0] - sliver, extent[1] + sliver]
                the_xticks = np.array([extent[0]])
            else:
                dx, = np.diff(extent[:2])/(x_amps_data.T.shape[1]-1)
                shifted_extent = [extent[0]-dx/2, extent[1]+dx/2]
                x_step = dx
                while 20 * x_step < np.diff(extent[:2]):
                    x_step = 2 * x_step
                the_xticks = np.arange(extent[0], extent[1]+dx, x_step)
            if np.diff(extent[2:]) <= 10**-6:
                shifted_extent.extend([extent[2] - sliver, extent[3] + sliver])
                the_yticks = np.array([extent[2]])
            else:
                dy, = np.diff(extent[2:])/(x_amps_data.T.shape[0]-1)
                shifted_extent.extend([extent[2]-dy/2, extent[3]+dy/2])
                y_step = dy
                while 20 * y_step < np.diff(extent[2:]):
                    y_step = 2 * y_step
                the_yticks = np.arange(extent[2], extent[3]+dy, y_step)
            im = axs.imshow(x_amps_data.T, interpolation="none", origin='lower',
                            cmap="viridis", extent=shifted_extent, aspect='auto')
            axs.set_xticks(the_xticks)
            axs.set_yticks(the_yticks)
            axs.set_ylabel("beta")
            axs.set_xlabel("alpha")
            tol_type = "rel" if use_relative else "abs"
            the_title = "A graph of the x wave amplitude\nfor the simple truncated reduced model,\n"
            the_title += f"with {x0=}, {y0=}, t_max={interval[1]}, grains={granularity}," + "\n"
            the_title += f"abs sim tol={atol}, rel sim tol={rtol}, {tol_type} analysis tol={tolerance}."
            fig.suptitle(the_title)
            fig.colorbar(im, ax=axs, label='x amplitude (0 if converges)')
            plt.show()
        if plot_y_amps:
            have_y_amps = ["y_amplitude" in item for item in metrics_list]
            y_amps = [metrics_list[index]["y_amplitude"]
                      if have_y_amps[index] else np.nan for index in metrics_range]
            y_amps = [0 if are_collapsed[index] else y_amps[index] for index in metrics_range]
            y_amps_data = griddata(coefficients_list, y_amps, desired_range, method="nearest")

            fig, axs = plt.subplots(1, 1, layout='constrained')
            if np.diff(extent[:2]) <= sliver:
                shifted_extent = [extent[0] - sliver, extent[1] + sliver]
                the_xticks = np.array([extent[0]])
            else:
                dx, = np.diff(extent[:2])/(y_amps_data.T.shape[1]-1)
                shifted_extent = [extent[0]-dx/2, extent[1]+dx/2]
                x_step = dx
                while 20 * x_step < np.diff(extent[:2]):
                    x_step = 2 * x_step
                the_xticks = np.arange(extent[0], extent[1]+dx, x_step)
            if np.diff(extent[2:]) <= 10**-6:
                shifted_extent.extend([extent[2] - sliver, extent[3] + sliver])
                the_yticks = np.array([extent[2]])
            else:
                dy, = np.diff(extent[2:])/(y_amps_data.T.shape[0]-1)
                shifted_extent.extend([extent[2]-dy/2, extent[3]+dy/2])
                y_step = dy
                while 20 * y_step < np.diff(extent[2:]):
                    y_step = 2 * y_step
                the_yticks = np.arange(extent[2], extent[3]+dy, y_step)
            im = axs.imshow(y_amps_data.T, interpolation="none", origin='lower',
                            cmap="viridis", extent=shifted_extent, aspect='auto')
            axs.set_xticks(the_xticks)
            axs.set_yticks(the_yticks)
            axs.set_ylabel("beta")
            axs.set_xlabel("alpha")
            tol_type = "rel" if use_relative else "abs"
            the_title = "A graph of the y wave amplitude\nfor the simple truncated reduced model,\n"
            the_title += f"with {x0=}, {y0=}, t_max={interval[1]}, grains={granularity}," + "\n"
            the_title += f"abs sim tol={atol}, rel sim tol={rtol}, {tol_type} analysis tol={tolerance}."
            fig.suptitle(the_title)
            fig.colorbar(im, ax=axs, label='y amplitude (0 if converges)')
            plt.show()
        if plot_collapse_moments:
            collapse_moments = [item["collapse_moment"]
                                if "collapse_moment" in item else np.nan for item in metrics_list]
            collapse_moments_data = griddata(coefficients_list, collapse_moments, desired_range, method="nearest")

            fig, axs = plt.subplots(1, 1, layout='constrained')
            if np.diff(extent[:2]) <= sliver:
                shifted_extent = [extent[0] - sliver, extent[1] + sliver]
                the_xticks = np.array([extent[0]])
            else:
                dx, = np.diff(extent[:2])/(collapse_moments_data.T.shape[1]-1)
                shifted_extent = [extent[0]-dx/2, extent[1]+dx/2]
                x_step = dx
                while 20 * x_step < np.diff(extent[:2]):
                    x_step = 2 * x_step
                the_xticks = np.arange(extent[0], extent[1]+dx, x_step)
            if np.diff(extent[2:]) <= 10**-6:
                shifted_extent.extend([extent[2] - sliver, extent[3] + sliver])
                the_yticks = np.array([extent[2]])
            else:
                dy, = np.diff(extent[2:])/(collapse_moments_data.T.shape[0]-1)
                shifted_extent.extend([extent[2]-dy/2, extent[3]+dy/2])
                y_step = dy
                while 20 * y_step < np.diff(extent[2:]):
                    y_step = 2 * y_step
                the_yticks = np.arange(extent[2], extent[3]+dy, y_step)
            im = axs.imshow(collapse_moments_data.T, interpolation="none", origin='lower',
                            cmap="viridis", extent=shifted_extent, aspect='auto')
            axs.set_xticks(the_xticks)
            axs.set_yticks(the_yticks)
            axs.set_ylabel("beta")
            axs.set_xlabel("alpha")
            tol_type = "rel" if use_relative else "abs"
            the_title = "A graph of \"the moment\" x and y converge\nfor the simple truncated reduced model,\n"
            the_title += f"with {x0=}, {y0=}, t_max={interval[1]}, grains={granularity}," + "\n"
            the_title += f"abs sim tol={atol}, rel sim tol={rtol}, {tol_type} analysis tol={tolerance}."
            fig.suptitle(the_title)
            fig.colorbar(im, ax=axs, label='moment of convergence within tolerances')
            plt.show()

    def f(t, u, coeffs):
        x, y = u
        alpha, beta = coeffs
        return np.array([
            x * x * y - x - beta * x + alpha,
            x - x * x * y
        ])

    ##################################
    base_initials = [0, 0]  # list of starting values of the variables; first part of the parameters.
    base_coefficients = [0, 0]  # list of coefficients for reaction speeds; second part of the parameters.
    interval = (0, 1000)  # cutoff point in time to stop the simulation at, or None for the default value of 50.
    granularity = 2000  # number of points in time to actually log the values at (not counting t=0),
    # or None to let the solver itself decide for us.

    vary_simultaneously = False  # whether to entrywise combine the variations (True) or Cartesian them (False)
    multiplicative = False  # whether to apply variations multiplicatively (True) or additively (False)
    variations_initials = [None, None]
    # variations_coefficients = [np.linspace(0, 1, 51)[1:], np.linspace(0, 1, 51)[1:]]
    variations_coefficients = [np.linspace(0, 800, 33)[1:] / 2000, np.linspace(0, 2000, 41)[1:] / 2000]

    ############################################

    tmp = np.array([int(multiplicative)])
    variations_initials = [[item] if isinstance(item, int | float) else item for item in variations_initials]
    variations_coefficients = [[item] if isinstance(item, int | float) else item for item in variations_coefficients]
    variations_initials = [tmp if item is None else np.array(item) for item in variations_initials]
    variations_coefficients = [tmp if item is None else np.array(item) for item in variations_coefficients]

    initials_length = len(base_initials)
    base_variables = base_initials + base_coefficients
    variations_variables = variations_initials + variations_coefficients

    if vary_simultaneously:
        parallel_length = max([len(item) for item in variations_variables])
        for index in range(len(variations_variables)):
            if len(variations_variables[index]) == 1:
                variations_variables[index] = np.repeat(variations_variables[index], parallel_length)
            elif len(variations_variables[index]) != parallel_length:
                raise ValueError("Arrays aren't all of the same length or singletons")
        variations_combined = np.transpose(np.array(variations_variables))
    else:
        variations_combined = cartesian_product(*variations_variables)

    tmp = np.array(base_variables)
    variables_list = [tmp * item if multiplicative else tmp + item for item in variations_combined]
    coefficients_list = [variables[initials_length:] for variables in variables_list]

    evaluations = None if granularity is None else np.linspace(interval[0], interval[1], granularity + 1)

    simulate_and_plot_metrics(f, base_initials, coefficients_list, interval, granularity, atol=10**-7, rtol=10**-6, tolerance=10**-4, evaluations=evaluations, use_relative=False,
                              plot_periods=True, plot_x_amps=True, plot_y_amps=True, plot_collapse_moments=True)

    # TODO: debug this; it's giving zero-only plots and stuff
