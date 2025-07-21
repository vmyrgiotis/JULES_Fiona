import numpy as np
from scipy.integrate import solve_ivp

from triffid.triffid import triffid_rhs, params as triffid_params
from rothc.equations import soil_carbon_rhs
from rothc.parameters import POOLS, C0_default

def run_coupled_model(t_span_weeks, initial_conditions):

    t0_weeks, tf_weeks = t_span_weeks
    dt = 1.0  
    n_steps = int((tf_weeks - t0_weeks) / dt) + 1
    time_weeks = np.arange(t0_weeks, tf_weeks + dt, dt)
    
    # Pre-allocate arrays
    triffid_out = np.zeros((4, n_steps))
    rothc_out = np.zeros((4, n_steps))
    
    triffid_out[:, 0] = initial_conditions['triffid']
    rothc_out[:, 0] = initial_conditions['rothc']
    
    base_drivers = {
        'T_soil': lambda t: 283.15,
        'moisture': lambda t: 0.5,
        's': lambda t: 0.5
    }
    
    for i in range(1, n_steps):
        t = time_weeks[i-1]
        
        triffid_state = triffid_out[:, i-1]
        triffid_out[:, i] = triffid_state + dt * triffid_rhs(t, triffid_state)
        
        nu = triffid_out[2:, i]  
        nu_mean = nu.mean()
        
        rothc_drivers = {
            **base_drivers,
            'nu': lambda t: nu_mean,
            'Lambda_c': lambda t: 0.001
        }
        
        rothc_state = rothc_out[:, i-1]
        rothc_out[:, i] = rothc_state + dt * soil_carbon_rhs(t, rothc_state, rothc_drivers)
    
    return {
        'time_weeks': time_weeks,
        'triffid': triffid_out,
        'rothc': rothc_out
    }