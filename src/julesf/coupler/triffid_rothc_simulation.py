"""
TRIFFID-RothC Coupled Simulation - Minimal Version
Fixed dimension and bounds issues
"""
import numpy as np

from julesf.triffid.run_triffid import triffid_rhs, params as triffid_params
from julesf.rothc.equations import soil_carbon_rhs
from julesf.rothc.run_rothc import seasonal_temp, seasonal_moisture

def solve_triffid_rothc_euler(t_span, triffid_init, rothc_init, dt=1.0):
    """
    Euler integration for TRIFFID-RothC coupling
    Fixed array dimensions and bounds checking
    """
    t0, tf = t_span
    n_steps = int((tf - t0) / dt) + 1
    t = np.linspace(t0, tf, n_steps)  # Use linspace for exact dimensions
    
    # Pre-allocate - match time array length exactly
    triffid_out = np.zeros((4, n_steps))
    rothc_out = np.zeros((4, n_steps))
    litterfall = np.zeros(n_steps)
    veg_cover = np.zeros(n_steps)
    
    # Initial conditions
    triffid_out[:, 0] = triffid_init
    rothc_out[:, 0] = rothc_init
    litterfall[0] = 0.0
    veg_cover[0] = triffid_init[2] + triffid_init[3]  # Initial cover
    
    for i in range(1, n_steps):
        t_curr = t[i-1]
        
        # TRIFFID step with bounds checking
        triffid_state = triffid_out[:, i-1]
        dtriffid_dt = triffid_rhs(t_curr, triffid_state)
        new_triffid = triffid_state + dt * dtriffid_dt
        
        # Apply bounds to prevent negative LAI (fixes power warning)
        new_triffid[0] = max(0.01, new_triffid[0])  # Lb_tree >= 0.01
        new_triffid[1] = max(0.01, new_triffid[1])  # Lb_grass >= 0.01
        new_triffid[2] = max(0, min(1, new_triffid[2]))  # nu_tree in [0,1]
        new_triffid[3] = max(0, min(1, new_triffid[3]))  # nu_grass in [0,1]
        
        triffid_out[:, i] = new_triffid
        
        # Extract coupling variables
        Lb_tree, Lb_grass = new_triffid[0], new_triffid[1]
        nu_tree, nu_grass = new_triffid[2], new_triffid[3]
        
        # Litterfall coupling
        L_tree = triffid_params['sigma1'][0] * Lb_tree
        L_grass = triffid_params['sigma1'][1] * Lb_grass
        Lambda_c = triffid_params['gamma_l'][0] * (L_tree + L_grass) * 0.01
        nu_total = max(0, min(1, nu_tree + nu_grass))
        
        litterfall[i] = Lambda_c
        veg_cover[i] = nu_total
        
        # RothC step with coupling
        rothc_drivers = {
            'T_soil': seasonal_temp(t_curr),
            's': seasonal_moisture(t_curr),
            'nu': nu_total,
            'Lambda_c': Lambda_c,
        }
        
        rothc_state = rothc_out[:, i-1]
        drothc_dt = soil_carbon_rhs(t_curr, rothc_state, rothc_drivers)
        rothc_out[:, i] = rothc_state + dt * drothc_dt
    
    return t, triffid_out, rothc_out, {'litterfall': litterfall, 'vegetation_cover': veg_cover}
