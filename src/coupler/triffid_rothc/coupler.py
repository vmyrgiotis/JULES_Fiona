import numpy as np
from scipy.integrate import solve_ivp

# Import from absolute paths to avoid confusion
from src.triffid.triffid import triffid_rhs, params as triffid_params
from src.rothc.equations import soil_carbon_rhs
from src.rothc.parameters import POOLS, C0_default

def run_coupled_model(t_span_weeks, initial_conditions):
    """
    Run TRIFFID-RothC coupled model with weekly timesteps
    
    Parameters
    ----------
    t_span_weeks : tuple (t0, tf)
        Start and end times in weeks
    initial_conditions : dict
        'triffid': [Lb1, Lb2, nu1, nu2]
        'rothc': [C_DPM, C_RPM, C_BIO, C_HUM]
    """
    # Convert weekly timespan to both years (TRIFFID) and days (RothC)
    t0_weeks, tf_weeks = t_span_weeks
    dt_weeks = 1.0  # 1 week timestep
    
    # Time conversions
    weeks_to_years = 1.0 / 52.0
    weeks_to_days = 7.0
    
    n_steps = int((tf_weeks - t0_weeks) / dt_weeks) + 1
    
    # Pre-allocate output arrays
    triffid_out = np.zeros((4, n_steps))  # 2 LAI + 2 nu
    rothc_out = np.zeros((4, n_steps))    # 4 carbon pools
    time_weeks = np.arange(t0_weeks, tf_weeks + dt_weeks, dt_weeks)
    
    # Set initial conditions
    triffid_out[:, 0] = initial_conditions['triffid']
    rothc_out[:, 0] = initial_conditions['rothc']
    
    # Setup baseline drivers for RothC
    base_drivers = {
        'T_soil': lambda t: 283.15,  # 10°C
        'moisture': lambda t: 0.5,    # 50% moisture
        's': lambda t: 0.5           # soil moisture factor
    }
    
    # Main timestepping loop
    for i in range(1, n_steps):
        # Current timestep in different units
        t_now_weeks = time_weeks[i-1]
        t_next_weeks = time_weeks[i]
        
        t_now_years = t_now_weeks * weeks_to_years
        t_next_years = t_next_weeks * weeks_to_years
        
        t_now_days = t_now_weeks * weeks_to_days
        t_next_days = t_next_weeks * weeks_to_days
        
        # 1. Run TRIFFID (in years)
        triffid_sol = solve_ivp(
            triffid_rhs,
            (t_now_years, t_next_years),
            triffid_out[:, i-1],
            method='RK45'
        )
        
        # 2. Get vegetation cover for RothC
        nu = triffid_sol.y[2:, -1]  # Last two states are nu values
        nu_mean = nu.mean()  # Use mean cover for RothC
        
        # 3. Update RothC drivers with new nu
        rothc_drivers = {
            **base_drivers,
            'nu': lambda t: nu_mean,
            'Lambda_c': lambda t: 0.001  # Constant litter input for now
        }
        
        # 4. Run RothC (in days)
        def rothc_wrapped(t, y):
            return soil_carbon_rhs(t, y, rothc_drivers)
        
        rothc_sol = solve_ivp(
            rothc_wrapped,
            (t_now_days, t_next_days),
            rothc_out[:, i-1],
            method='RK45'
        )
        
        # Store results
        triffid_out[:, i] = triffid_sol.y[:, -1]
        rothc_out[:, i] = rothc_sol.y[:, -1]
    
    return {
        'time_weeks': time_weeks,
        'triffid': triffid_out,
        'rothc': rothc_out
    }