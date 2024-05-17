import scipy.integrate as scint

def simulate(f, initials_list, coefficients_list, interval, evaluations = None, atol=10**-7, rtol=10**-6):
    solution_list = list()

    for iteration in range(len(initials_list)):
        sol = scint.solve_ivp(fun=f,
                              t_span=interval,
                              y0=initials_list[iteration],
                              method="Radau",
                              t_eval=evaluations,
                              args=(coefficients_list[iteration],),
                              atol=atol,
                              rtol=rtol)
        
        solution_list.append(sol)
    
    return solution_list
