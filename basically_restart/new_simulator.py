import scipy.integrate as scint

def simulate(f, initials_list, coefficients_list, interval_list, evaluations_list = None):
    solution_list = list()



    for iteration in range(len(initials_list)):
        t_eval = None if evaluations_list is None else evaluations_list[iteration]

        sol = scint.solve_ivp(fun=f,
                              t_span=interval_list[iteration],
                              y0=initials_list[iteration],
                              method="Radau",
                              t_eval=t_eval,
                              args=(coefficients_list[iteration],))
        
        solution_list.append(sol)
    
    return (f, initials_list, coefficients_list, interval_list, evaluations_list, solution_list)
