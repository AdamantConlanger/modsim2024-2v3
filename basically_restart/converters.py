import math


def simple_unreduced_to_simple_reduced(initials, coefficients, interval=None, evaluations=None, magnitudes_list=None):
    p0_old, q0_old, x0_old, y0_old, r0_old, s0_old = initials
    k1_old, k2_old, k3_old, k4_old = coefficients

    factors = (math.sqrt(k3_old / (k2_old * q0_old)), k2_old * q0_old)

    x0 = x0_old * factors[0]
    y0 = y0_old * factors[0]
    k1 = k1_old * p0_old * factors[0] / factors[1]
    k2 = k4_old / factors[1]

    initials_new = [x0, y0]
    coefficients_new = [k1, k2]

    if interval is None:
        interval_new = None
    else:
        tlower_old, tupper_old = interval
        tupper = (tupper_old - tlower_old) * factors[1]
        interval_new = (0, tupper)

    if evaluations is None:
        evaluations_new = None
    else:
        evaluations_new = [(evaluation_old - tlower_old) * factors[1] for evaluation_old in evaluations]

    if magnitudes_list is None:
        magnitudes_new_list = None
    else:
        magnitudes_new_list = list()
        for magnitudes_old in magnitudes_list:
            p_old, q_old, x_old, y_old, r_old, s_old = magnitudes_old
            x = x_old * factors[0]
            y = y_old * factors[0]
            magnitudes_new_list.append([x, y])

    return (initials_new, coefficients_new, interval_new, evaluations_new, magnitudes_new_list, initials, coefficients, interval, evaluations, magnitudes_list)


def simple_unreduced_to_simple_reduced_list(initials_list, coefficients_list, interval_list=None, evaluations_list=None, magnitudes_list_list=None):
    initials_new_list = list()
    coefficients_new_list = list()
    interval_new_list = None if interval_list is None else list()
    evaluations_new_list = None if evaluations_list is None else list()
    magnitudes_new_list_list = None if magnitudes_list_list is None else list()
    for iteration in range(len(initials_list)):
        interval = None if interval_list is None else interval_list[iteration]
        evaluations = None if evaluations_list is None else evaluations_list[iteration]
        magnitudes_list = None if magnitudes_list_list is None else magnitudes_list_list[iteration]
        result = simple_unreduced_to_simple_reduced(
            initials_list[iteration], coefficients_list[iteration], interval, evaluations, magnitudes_list)
        initials_new_list.append(result[0])
        coefficients_new_list.append(result[1])
        if interval_list is not None:
            interval_new_list.append(result[2])
        if evaluations_list is not None:
            evaluations_new_list.append(result[3])
        if magnitudes_list_list is not None:
            magnitudes_new_list_list.append(result[4])
    return (initials_new_list, coefficients_new_list, interval_new_list, evaluations_new_list, magnitudes_new_list_list, initials_list, coefficients_list, interval_list, evaluations_list, magnitudes_list_list)


def simple_reduced_to_simple_unreduced(initials_unred, coefficients_unred, interval_unred, initials, coefficients, interval=None, evaluations=None, magnitudes_list=None):
    p0_unred, q0_unred, x0_unred, y0_unred, r0_unred, s0_unred = initials_unred
    k1_unred, k2_unred, k3_unred, k4_unred = coefficients_unred
    x0_old, y0_old = initials
    k1_old, k2_old = coefficients

    factors = (math.sqrt(k2_unred * q0_unred / k3_unred), k2_unred * q0_unred)

    initials_new = initials_unred
    coefficients_new = coefficients_unred
    p0, q0, x0, y0, r0, s0 = initials_new
    k1, k2, k3, k4 = coefficients_new

    if interval is None:
        interval_new = None
        tlower = 0
    elif interval_unred is None:
        # raise ValueError("Can't convert intervals because unreduced interval is None.")
        print("Can't convert intervals because unreduced interval is None; assuming tlower = 0")
        tlower_old, tupper_old = interval
        tlower = 0
        tupper = (tupper_old - tlower_old) / factors[1]
        interval_new = (tlower, tupper)
    else:
        tlower_old, tupper_old = interval
        tlower_unred, tupper_unred = interval_unred
        interval_new = interval_unred
        tlower, tupper = interval_new

    if evaluations is None:
        evaluations_new = None
    else:
        evaluations_new = [tlower + (evaluation_old - tlower_old) / factors[1] for evaluation_old in evaluations]

    if magnitudes_list is None:
        magnitudes_new_list = None
    else:
        magnitudes_new_list = list()
        for magnitudes_old in magnitudes_list:
            x_old, y_old = magnitudes_old
            p = p0
            q = q0
            x = x_old * factors[0]
            y = y_old * factors[0]
            r = None
            s = None
            magnitudes_new_list.append([p, q, x, y, r, s])

    return (initials_new, coefficients_new, interval_new, evaluations_new, magnitudes_new_list, initials, coefficients, interval, evaluations, magnitudes_list)


def simple_reduced_to_simple_unreduced_list(initials_unred_list, coefficients_unred_list, interval_unred_list, initials_list, coefficients_list, interval_list=None, evaluations_list=None, magnitudes_list_list=None):
    initials_new_list = list()
    coefficients_new_list = list()
    interval_new_list = None if interval_list is None else list()
    evaluations_new_list = None if evaluations_list is None else list()
    magnitudes_new_list_list = None if magnitudes_list_list is None else list()
    for iteration in range(len(initials_list)):
        initials_unred = initials_unred_list[iteration]
        coefficients_unred = coefficients_unred_list[iteration]
        interval_unred = interval_unred_list[iteration]
        interval = None if interval_list is None else interval_list[iteration]
        evaluations = None if evaluations_list is None else evaluations_list[iteration]
        magnitudes_list = None if magnitudes_list_list is None else magnitudes_list_list[iteration]
        result = simple_reduced_to_simple_unreduced(initials_unred, coefficients_unred, interval_unred,
            initials_list[iteration], coefficients_list[iteration], interval, evaluations, magnitudes_list)
        initials_new_list.append(result[0])
        coefficients_new_list.append(result[1])
        if interval_list is not None:
            interval_new_list.append(result[2])
        if evaluations_list is not None:
            evaluations_new_list.append(result[3])
        if magnitudes_list_list is not None:
            magnitudes_new_list_list.append(result[4])
    return (initials_new_list, coefficients_new_list, interval_new_list, evaluations_new_list, magnitudes_new_list_list, initials_list, coefficients_list, interval_list, evaluations_list, magnitudes_list_list)

    

