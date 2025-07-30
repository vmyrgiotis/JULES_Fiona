#simulation.py

from scipy.integrate import solve_ivp
import numpy as np
from .equations import triffid_rhs

def solve_triffid(t_span, initial_conditions, params, npp_external=None, method='RK45'):
    """
    Solve TRIFFID equations
    
    Returns:
        t: time array
        Lb: balanced LAI [n_pfts, n_times]  
        nu: fractional cover [n_pfts, n_times]
    """
    y0 = np.concatenate([initial_conditions['Lb_init'], 
                         initial_conditions['nu_init']])
    
    t_eval = np.linspace(t_span[0], t_span[1], 501)
    
    sol = solve_ivp(
        lambda t, y: triffid_rhs(t, y, params, npp_external),
        t_span, y0, t_eval=t_eval, method=method
    )
    
    return sol.t, sol.y[:2], sol.y[2:]