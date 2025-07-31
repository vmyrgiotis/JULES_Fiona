"""
Wrapper for slow models (TRIFFID + RothC) at weekly timestep
"""
import numpy as np

from julesf.triffid.equations import triffid_rhs
from julesf.triffid.parameters import PFT_PARAMS
from julesf.rothc.equations import soil_carbon_rhs
from julesf.coupler.coupling_utils import calculate_litterfall, weekly_mean_soil_moisture

def run_slow_models_week(fast_results, triffid_state, rothc_state, week_num=0):
    """
    Run TRIFFID + RothC for 1 week timestep
    
    Args:
        fast_results: Output from run_fast_models_week()
        triffid_state: [Lb_tree, Lb_grass, nu_tree, nu_grass]
        rothc_state: [C_dpm, C_rpm, C_bio, C_hum]
        week_num: Week number for seasonal forcing
        
    Returns:
        Dict with updated states and coupling variables
    """
    # 1. Run TRIFFID with real NPP from physiology
    npp_week = fast_results['NPP_total']
    triffid_new = run_triffid_with_npp(npp_week, triffid_state, dt_weeks=1)
    
    # 2. Calculate litterfall from new TRIFFID state
    litterfall = calculate_litterfall(triffid_new, PFT_PARAMS)
    
    # 3. Extract vegetation cover for RothC
    nu_total = triffid_new[2] + triffid_new[3]
    
    # 4. Get mean soil moisture for RothC
    s = weekly_mean_soil_moisture(fast_results['soil_theta'])
    
    # 5. Get mean soil temperature
    T_soil_mean = np.mean(fast_results['soil_T'][0, :]) - 273.15  # Convert K to °C
    
    # 6. Run RothC
    rothc_drivers = {
        'T_soil': T_soil_mean,
        's': s,
        'nu': nu_total,
        'Lambda_c': litterfall,
    }
    
    t = week_num * 7  # Current time in days
    drothc_dt = soil_carbon_rhs(t, rothc_state, rothc_drivers)
    rothc_new = rothc_state + 7 * drothc_dt  # 7 days timestep
    
    return {
        'triffid_new': triffid_new,
        'rothc_new': rothc_new,
        'coupling_vars': {
            'litterfall': litterfall,
            'nu_total': nu_total
        }
    }

def run_triffid_with_npp(npp_kg_c_week, triffid_state, dt_weeks=1):
    """Run TRIFFID with real NPP from physiology"""
    dt_years = dt_weeks / 52  # Convert weeks to years for TRIFFID
    
    # Create parameters for TRIFFID
    params_combined = {
        **PFT_PARAMS, 
        'competition': {'c': np.array([[1.0, 1.0], [0.0, 1.0]])}
    }
    
    # Convert weekly NPP to yearly for TRIFFID
    npp_yearly = npp_kg_c_week * 52
    npp_external = np.array([npp_yearly, npp_yearly])  # Same for all PFTs
    
    # Calculate derivatives
    dtriffid_dt = triffid_rhs(0, triffid_state, params_combined, npp_external=npp_external)
    
    # Euler step (dt in years)
    triffid_new = triffid_state + dt_years * dtriffid_dt
    
    # Apply bounds for stability
    triffid_new[0] = max(0.01, triffid_new[0])  # Lb_tree
    triffid_new[1] = max(0.01, triffid_new[1])  # Lb_grass
    triffid_new[2] = max(0, min(1, triffid_new[2]))  # nu_tree
    triffid_new[3] = max(0, min(1, triffid_new[3]))  # nu_grass
    
    return triffid_new