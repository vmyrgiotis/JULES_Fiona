# simulation.py

import numpy as np
from scipy.integrate import solve_ivp

try: 
    from equations import soil_carbon_rhs
except ModuleNotFoundError:
    from .equations import soil_carbon_rhs

def run_soil_model(t_span, C0, drivers_ts,
                   texture="loam", temp_fun="Q10",
                   pft="tree",
                   method="RK45", dt_out=1.0):
    """
    Integrate the soil-carbon model.

    Parameters
    ----------
    t_span     : tuple (t0, tf) in days
    C0         : array_like, initial [C_DPM, C_RPM, C_BIO, C_HUM]
    drivers_ts : dict of callables:
                 'Lambda_c', 'T_soil', 's', 'nu' → f(t) → scalar
    texture    : soil texture string for β_R lookup
    temp_fun   : "Q10" or "RothC"
    dt_out     : output interval (days)

    Returns
    -------
    t  : 1d array of times
    C  : 2d array, shape (4, len(t)), pool values over time
    """
    def rhs(t, C):
        drv = {k: f(t) for k, f in drivers_ts.items()}
        return soil_carbon_rhs(t, C, drv,
                               texture=texture,
                               temp_fun=temp_fun,
                               pft=pft)

    t_eval = np.arange(t_span[0], t_span[1] + dt_out, dt_out)
    sol = solve_ivp(rhs, t_span, C0,
                    method=method,
                    t_eval=t_eval)
    return sol.t, sol.y
