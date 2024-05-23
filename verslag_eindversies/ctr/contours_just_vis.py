def perform_program(simulate, cartesian_product):
    import scipy.integrate as scint
    from scipy.interpolate import griddata
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    import json

    plt.rcParams.update({'font.size': 18})

    def plot_metrics(filename, f, initials, coefficients_list, interval, granularity, atol=10**-7, rtol=10**-6, tolerance=10**-4, evaluations=None, use_relative=False, plot_periods=True, plot_x_amps=False, plot_y_amps=False, plot_collapse_moments=False, autolevels=False):
        x0, y0 = initials
        alpha_max = max([item[0] for item in coefficients_list])
        alpha_min = min([item[0] for item in coefficients_list])
        beta_max = max([item[1] for item in coefficients_list])
        beta_min = min([item[1] for item in coefficients_list])
        alpha_len = len(set([item[0] for item in coefficients_list]))
        beta_len = len(set([item[1] for item in coefficients_list]))
        extent = (alpha_min, alpha_max, beta_min, beta_max)

        with open(filename, "r") as the_file:
            metrics_list = json.load(the_file)

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

            periods_maximum = max(periods)
            nonzero_periods = np.array(periods)[np.array(periods) >= sliver]
            if len(nonzero_periods) == 0:
                print("couldn't create period contours; all periods are the same.")
            elif min(nonzero_periods) + sliver >= periods_maximum:
                print("couldn't create period contours; all periods are the same.")
            else:
                periods_minimum = min(nonzero_periods)
                num_of_contours = 20
                if autolevels:
                    levels = np.exp(np.linspace(math.log(periods_minimum), math.log(periods_maximum), num_of_contours))
                else:
                    levels = num_of_contours

                fig, axs = plt.subplots(1, 1, layout='constrained')
                if np.diff(extent[:2]) <= sliver:
                    shifted_extent = [extent[0] - sliver, extent[1] + sliver]
                    the_xticks = np.array([extent[0]])
                else:
                    dx, = np.diff(extent[:2])/(periods_data.T.shape[1]-1)
                    shifted_extent = [extent[0]-dx/2, extent[1]+dx/2]
                    x_step = dx
                    while 10 * x_step < np.diff(extent[:2]):
                        x_step = 2 * x_step
                    the_xticks = np.arange(extent[0], extent[1]+dx, x_step)
                if np.diff(extent[2:]) <= 10**-6:
                    shifted_extent.extend([extent[2] - sliver, extent[3] + sliver])
                    the_yticks = np.array([extent[2]])
                else:
                    dy, = np.diff(extent[2:])/(periods_data.T.shape[0]-1)
                    shifted_extent.extend([extent[2]-dy/2, extent[3]+dy/2])
                    y_step = dy
                    while 10 * y_step < np.diff(extent[2:]):
                        y_step = 2 * y_step
                    the_yticks = np.arange(extent[2], extent[3]+dy, y_step)
                contour_fills = axs.contourf(periods_data.T, levels, origin='lower',
                                             cmap="viridis", extent=shifted_extent)
                contours = axs.contour(periods_data.T, contour_fills.levels,
                                       origin='lower', colors='k', extent=shifted_extent)
                axs.set_xticks(the_xticks)
                axs.set_yticks(the_yticks)
                axs.set_ylabel("beta*")
                axs.set_xlabel("alpha*")
                tol_type = "rel" if use_relative else "abs"
                the_title = "Contourplot van de oscillatieperiode\nvoor de gereduceerde uitbreiding,\n"
                the_title += f"x*0={x0}, y*0={y0}, t_max={interval[1]}, grains={granularity}," + "\n"
                the_title += f"abs sim tol={atol}, rel sim tol={rtol}, {tol_type} analysetol={tolerance}."
                fig.suptitle(the_title)
                fig.colorbar(contour_fills, ax=axs, label='oscillatieperiode (0 bij convergentie)')
                plt.show()
        if plot_x_amps:
            have_x_amps = ["x_amplitude" in item for item in metrics_list]
            x_amps = [metrics_list[index]["x_amplitude"]
                      if have_x_amps[index] else np.nan for index in metrics_range]
            x_amps = [0 if are_collapsed[index] else x_amps[index] for index in metrics_range]
            x_amps_data = griddata(coefficients_list, x_amps, desired_range, method="nearest")

            x_amps_maximum = max(x_amps)
            nonzero_x_amps = np.array(x_amps)[np.array(x_amps) >= sliver]
            if len(nonzero_x_amps) == 0:
                print("couldn't create x amplitude contours; all x amplitudes are the same.")
            elif min(nonzero_x_amps) + sliver >= x_amps_maximum:
                print("couldn't create x amplitude contours; all x amplitudes are the same.")
            else:
                x_amps_minimum = min(nonzero_x_amps)
                num_of_contours = 20
                if autolevels:
                    levels = np.exp(np.linspace(math.log(x_amps_minimum), math.log(x_amps_maximum), num_of_contours))
                else:
                    levels = num_of_contours

                fig, axs = plt.subplots(1, 1, layout='constrained')
                if np.diff(extent[:2]) <= sliver:
                    shifted_extent = [extent[0] - sliver, extent[1] + sliver]
                    the_xticks = np.array([extent[0]])
                else:
                    dx, = np.diff(extent[:2])/(x_amps_data.T.shape[1]-1)
                    shifted_extent = [extent[0]-dx/2, extent[1]+dx/2]
                    x_step = dx
                    while 10 * x_step < np.diff(extent[:2]):
                        x_step = 2 * x_step
                    the_xticks = np.arange(extent[0], extent[1]+dx, x_step)
                if np.diff(extent[2:]) <= 10**-6:
                    shifted_extent.extend([extent[2] - sliver, extent[3] + sliver])
                    the_yticks = np.array([extent[2]])
                else:
                    dy, = np.diff(extent[2:])/(x_amps_data.T.shape[0]-1)
                    shifted_extent.extend([extent[2]-dy/2, extent[3]+dy/2])
                    y_step = dy
                    while 10 * y_step < np.diff(extent[2:]):
                        y_step = 2 * y_step
                    the_yticks = np.arange(extent[2], extent[3]+dy, y_step)
                contour_fills = axs.contourf(x_amps_data.T, levels, origin='lower',
                                             cmap="viridis", extent=shifted_extent)
                contours = axs.contour(x_amps_data.T, contour_fills.levels,
                                       origin='lower', colors='k', extent=shifted_extent)
                axs.set_xticks(the_xticks)
                axs.set_yticks(the_yticks)
                axs.set_ylabel("beta*")
                axs.set_xlabel("alpha*")
                tol_type = "rel" if use_relative else "abs"
                the_title = "Contourplot van de amplitude van golven in x*\nvoor de gereduceerde uitbreiding,\n"
                the_title += f"x*0={x0}, y*0={y0}, t_max={interval[1]}, grains={granularity}," + "\n"
                the_title += f"abs sim tol={atol}, rel sim tol={rtol}, {tol_type} analysetol={tolerance}."
                fig.suptitle(the_title)
                fig.colorbar(contour_fills, ax=axs, label='x*-amplitude (0 bij convergentie)')
                plt.show()
        if plot_y_amps:
            have_y_amps = ["y_amplitude" in item for item in metrics_list]
            y_amps = [metrics_list[index]["y_amplitude"]
                      if have_y_amps[index] else np.nan for index in metrics_range]
            y_amps = [0 if are_collapsed[index] else y_amps[index] for index in metrics_range]
            y_amps_data = griddata(coefficients_list, y_amps, desired_range, method="nearest")

            y_amps_maximum = max(y_amps)
            nonzero_y_amps = np.array(y_amps)[np.array(y_amps) >= sliver]
            if len(nonzero_y_amps) == 0:
                print("couldn't create y amplitude contours; all y amplitudes are the same.")
            elif min(nonzero_y_amps) + sliver >= y_amps_maximum:
                print("couldn't create y amplitude contours; all y amplitudes are the same.")
            else:
                y_amps_minimum = min(nonzero_y_amps)
                num_of_contours = 20
                if autolevels:
                    levels = np.exp(np.linspace(math.log(y_amps_minimum), math.log(y_amps_maximum), num_of_contours))
                else:
                    levels = num_of_contours

                fig, axs = plt.subplots(1, 1, layout='constrained')
                if np.diff(extent[:2]) <= sliver:
                    shifted_extent = [extent[0] - sliver, extent[1] + sliver]
                    the_xticks = np.array([extent[0]])
                else:
                    dx, = np.diff(extent[:2])/(y_amps_data.T.shape[1]-1)
                    shifted_extent = [extent[0]-dx/2, extent[1]+dx/2]
                    x_step = dx
                    while 10 * x_step < np.diff(extent[:2]):
                        x_step = 2 * x_step
                    the_xticks = np.arange(extent[0], extent[1]+dx, x_step)
                if np.diff(extent[2:]) <= 10**-6:
                    shifted_extent.extend([extent[2] - sliver, extent[3] + sliver])
                    the_yticks = np.array([extent[2]])
                else:
                    dy, = np.diff(extent[2:])/(y_amps_data.T.shape[0]-1)
                    shifted_extent.extend([extent[2]-dy/2, extent[3]+dy/2])
                    y_step = dy
                    while 10 * y_step < np.diff(extent[2:]):
                        y_step = 2 * y_step
                    the_yticks = np.arange(extent[2], extent[3]+dy, y_step)
                contour_fills = axs.contourf(y_amps_data.T, levels, origin='lower',
                                             cmap="viridis", extent=shifted_extent)
                contours = axs.contour(y_amps_data.T, contour_fills.levels,
                                       origin='lower', colors='k', extent=shifted_extent)
                axs.set_xticks(the_xticks)
                axs.set_yticks(the_yticks)
                axs.set_ylabel("beta*")
                axs.set_xlabel("alpha*")
                tol_type = "rel" if use_relative else "abs"
                the_title = "Contourplot van de amplitude van golven in y*\nvoor de gereduceerde uitbreiding,\n"
                the_title += f"x*0={x0}, y*0={y0}, t_max={interval[1]}, grains={granularity}," + "\n"
                the_title += f"abs sim tol={atol}, rel sim tol={rtol}, {tol_type} analysetol={tolerance}."
                fig.suptitle(the_title)
                fig.colorbar(contour_fills, ax=axs, label='y*-amplitude (0 bij convergentie)')
                plt.show()
        if plot_collapse_moments:
            collapse_moments = [item["collapse_moment"]
                                if "collapse_moment" in item else np.nan for item in metrics_list]
            collapse_moments_data = griddata(coefficients_list, collapse_moments, desired_range, method="nearest")

            collapse_moments_maximum = max(collapse_moments)
            nonzero_collapse_moments = np.array(collapse_moments)[np.array(collapse_moments) >= sliver]
            if len(nonzero_collapse_moments) == 0:
                print("couldn't create collapse moment contours; all collapse moments are the same.")
            elif min(nonzero_collapse_moments) + sliver >= collapse_moments_maximum:
                print("couldn't create collapse moment contours; all collapse moments are the same.")
            else:
                collapse_moments_minimum = min(nonzero_collapse_moments)
                num_of_contours = 20
                if autolevels:
                    the_log_max = math.log(collapse_moments_maximum)
                    the_log_min = math.log(collapse_moments_minimum)
                    levels = np.exp(np.linspace(the_log_min, the_log_max, num_of_contours))
                else:
                    levels = num_of_contours

                fig, axs = plt.subplots(1, 1, layout='constrained')
                if np.diff(extent[:2]) <= sliver:
                    shifted_extent = [extent[0] - sliver, extent[1] + sliver]
                    the_xticks = np.array([extent[0]])
                else:
                    dx, = np.diff(extent[:2])/(collapse_moments_data.T.shape[1]-1)
                    shifted_extent = [extent[0]-dx/2, extent[1]+dx/2]
                    x_step = dx
                    while 10 * x_step < np.diff(extent[:2]):
                        x_step = 2 * x_step
                    the_xticks = np.arange(extent[0], extent[1]+dx, x_step)
                if np.diff(extent[2:]) <= 10**-6:
                    shifted_extent.extend([extent[2] - sliver, extent[3] + sliver])
                    the_yticks = np.array([extent[2]])
                else:
                    dy, = np.diff(extent[2:])/(collapse_moments_data.T.shape[0]-1)
                    shifted_extent.extend([extent[2]-dy/2, extent[3]+dy/2])
                    y_step = dy
                    while 10 * y_step < np.diff(extent[2:]):
                        y_step = 2 * y_step
                    the_yticks = np.arange(extent[2], extent[3]+dy, y_step)
                contour_fills = axs.contourf(collapse_moments_data.T, levels, origin='lower',
                                             cmap="viridis", extent=shifted_extent)
                contours = axs.contour(collapse_moments_data.T, contour_fills.levels,
                                       origin='lower', colors='k', extent=shifted_extent)
                axs.set_xticks(the_xticks)
                axs.set_yticks(the_yticks)
                axs.set_ylabel("beta*")
                axs.set_xlabel("alpha*")
                tol_type = "rel" if use_relative else "abs"
                the_title = "Contourplot van afschatting van convergentiemoment van x* en y*\n"
                the_title += "voor de gereduceerde uitbreiding,\n"
                the_title += f"x*0={x0}, y*0={y0}, t_max={interval[1]}, grains={granularity}," + "\n"
                the_title += f"abs sim tol={atol}, rel sim tol={rtol}, {tol_type} analysetol={tolerance}."
                fig.suptitle(the_title)
                fig.colorbar(contour_fills, ax=axs, label='convergentiemoment binnen tolerantie')
                plt.show()

    def f(t, u, coeffs):
        x, y = u
        alpha, beta = coeffs
        return np.array([
            x * x * x * y - x - beta * x + alpha,
            x - x * x * x * y
        ])

    ##################################
    base_initials = [0, 0]  # list of starting values of the variables; first part of the parameters.
    base_coefficients = [0, 0]  # list of coefficients for reaction speeds; second part of the parameters.
    interval = (0, 5000)  # cutoff point in time to stop the simulation at, or None for the default value of 50.
    granularity = 5000  # number of points in time to actually log the values at (not counting t=0),
    # or None to let the solver itself decide for us.

    vary_simultaneously = False  # whether to entrywise combine the variations (True) or Cartesian them (False)
    multiplicative = False  # whether to apply variations multiplicatively (True) or additively (False)
    variations_initials = [None, None]
    # variations_coefficients = [np.linspace(0, 1, 51)[1:], np.linspace(0, 1, 51)[1:]]
    variations_coefficients = [np.linspace(0, 250, 51)[1:] / 200, np.linspace(0, 400, 41)[1:] / 200]

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

    plot_metrics("ctr_metrics.json", f, base_initials, coefficients_list, interval, granularity, atol=10**-7, rtol=10**-6, tolerance=10**-3, evaluations=evaluations, use_relative=False,
                 plot_periods=True, plot_x_amps=True, plot_y_amps=True, plot_collapse_moments=True, autolevels=False)
