import numpy as np

# Updated imports for new TRIFFID structure
from julesf.triffid.equations import triffid_rhs
from julesf.triffid.parameters import PFT_PARAMS
from julesf.rothc.equations import soil_carbon_rhs
from julesf.rothc.run_rothc import seasonal_temp, seasonal_moisture

def solve_triffid_rothc_euler(t_span, triffid_init, rothc_init, dt=1.0):
    """
    Euler integration for TRIFFID-RothC coupling
    """
    t0, tf = t_span
    n_steps = int((tf - t0) / dt) + 1
    t = np.linspace(t0, tf, n_steps)
    
    # Pre-allocate
    triffid_out = np.zeros((4, n_steps))
    rothc_out = np.zeros((4, n_steps))
    litterfall = np.zeros(n_steps)
    veg_cover = np.zeros(n_steps)
    
    # Initial conditions
    triffid_out[:, 0] = triffid_init
    rothc_out[:, 0] = rothc_init
    litterfall[0] = 0.0
    veg_cover[0] = triffid_init[2] + triffid_init[3]
    
    # Prepare TRIFFID parameters in the format expected by triffid_rhs
    params_combined = {**PFT_PARAMS, 'competition': {'c': np.array([[1.0, 1.0], [0.0, 1.0]])}}
    
    for i in range(1, n_steps):
        t_curr = t[i-1]
        
        # TRIFFID step 
        triffid_state = triffid_out[:, i-1]
        dtriffid_dt = triffid_rhs(t_curr, triffid_state, params_combined, npp_external=None)
        new_triffid = triffid_state + dt * dtriffid_dt
        
        # Apply bounds
        new_triffid[0] = max(0.01, new_triffid[0])
        new_triffid[1] = max(0.01, new_triffid[1])
        new_triffid[2] = max(0, min(1, new_triffid[2]))
        new_triffid[3] = max(0, min(1, new_triffid[3]))
        
        triffid_out[:, i] = new_triffid
        
        # Extract coupling variables 
        Lb_tree, Lb_grass = new_triffid[0], new_triffid[1]
        nu_tree, nu_grass = new_triffid[2], new_triffid[3]
        
        # Litterfall coupling 
        sigma1_tree = PFT_PARAMS['broadleaf_tree']['sigma1']
        sigma1_grass = PFT_PARAMS['C3_grass']['sigma1']
        gamma_l_tree = PFT_PARAMS['broadleaf_tree']['gamma_l']
        gamma_l_grass = PFT_PARAMS['C3_grass']['gamma_l']
        
        L_tree = sigma1_tree * Lb_tree
        L_grass = sigma1_grass * Lb_grass
        Lambda_c = gamma_l_tree * L_tree + gamma_l_grass * L_grass  # kg C m⁻² yr⁻¹
        Lambda_c = Lambda_c / 365.25  # Convert to daily rate for coupling
        
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
