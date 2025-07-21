import numpy as np
from scipy.integrate import solve_ivp

from triffid.triffid import triffid_rhs, params as triffid_params
from rothc.equations import soil_carbon_rhs 
from rothc.parameters import POOLS, C0_default

def run_coupled_model(t_span, initial_conditions, drivers_ts, dt=1.0):
    """
    Run coupled TRIFFID-RothC simulation
    
    Parameters
    ----------
    t_span : tuple (t0, tf) 
        Start and end time in days
    initial_conditions : dict
        'triffid': [Lb1, Lb2, nu1, nu2]  # LAI and cover fractions
        'rothc': [C_DPM, C_RPM, C_BIO, C_HUM]  # Carbon pools
    drivers_ts : dict 
        'T_soil': callable t -> temp in K
        'moisture': callable t -> soil moisture (0-1)
    dt : float
        Coupling timestep in days
    """
    # initial states
    triffid_state = initial_conditions['triffid']
    rothc_state = initial_conditions['rothc']
    
    # time points
    t = np.arange(t_span[0], t_span[1]+dt, dt)
    n_steps = len(t)
    
    # pre-allocate output arrays
    triffid_out = np.zeros((4, n_steps))  # 4 = 2 LAI + 2 nu
    rothc_out = np.zeros((4, n_steps))    # 4 carbon pools
    
    # store initial conditions
    triffid_out[:,0] = triffid_state
    rothc_out[:,0] = rothc_state
    
    # main loop
    for i in range(1, n_steps):
        t_now = t[i-1]
        t_next = t[i]
        step_span = (t_now, t_next)
        
        # 1. Run TRIFFID forward one step
        triffid_sol = solve_ivp(
            triffid_rhs,
            step_span,
            triffid_out[:,i-1],
            method='RK45'
        )
        
        # 2. Extract vegetation cover (nu) for RothC
        nu_new = triffid_sol.y[2:,-1]  # Last two states are nu
        
        # 3. Update RothC drivers with new nu
        rothc_drivers = {
            **drivers_ts,
            'nu': lambda t: nu_new.mean(),  # Mean cover fraction
            'Lambda_c': lambda t: 0.001     # Constant litter input for now
        }
        
        # 4. Run RothC forward one step
        def rothc_wrapped(t, y):
            drivers = {k: v(t) for k,v in rothc_drivers.items()}
            return soil_carbon_rhs(t, y, drivers)
            
        rothc_sol = solve_ivp(
            rothc_wrapped, 
            step_span,
            rothc_out[:,i-1],
            method='RK45'
        )
        
        # 5. Store results
        triffid_out[:,i] = triffid_sol.y[:,-1]
        rothc_out[:,i] = rothc_sol.y[:,-1]
    
    return {
        'time': t,
        'triffid': triffid_out,
        'rothc': rothc_out
    }