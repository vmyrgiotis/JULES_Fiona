# simulation.py

import numpy as np
from scipy.integrate import solve_ivp
from equations import ebm_rhs


def run_ebm(t_span, T0, drivers, method="RK45", dt_out=0.1):
    def rhs(t, T):
        return ebm_rhs(t, T, drivers)

    t_eval = np.arange(t_span[0], t_span[1]+dt_out, dt_out)
    sol = solve_ivp(rhs, t_span, [T0], method=method, t_eval=t_eval)
    return sol.t, sol.y[0]
