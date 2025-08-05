import numpy as np
import time

from julesf.coupler.fast_models_wrapper import run_fast_models_week
from julesf.coupler.slow_models_wrapper import run_slow_models_week
from julesf.coupler.coupling_utils import extract_vegetation_vars
from julesf.rothc.parameters import POOLS, C0_default
from julesf.soil.parameters import SOIL_LAYERS, INITIAL_CONDITIONS

def jules_master_coupler(weeks=52, triffid_init=None, rothc_init=None, soil_init=None, 
                        external_drivers=None):
    """
    JULES master coupler with optional external meteorological drivers
    
    Args:
        weeks: Number of weeks to simulate
        triffid_init: Initial TRIFFID state [Lb_tree, Lb_grass, nu_tree, nu_grass]
        rothc_init: Initial RothC state [C_dpm, C_rpm, C_bio, C_hum]
        soil_init: Initial soil state {'theta': array, 'T_soil': array}
        external_drivers: Optional external forcing in format:
                         {
                             'era5_forcing': {
                                 'Tair': lambda t: temp_at_time_t,
                                 'Sdn': lambda t: solar_at_time_t,
                                 'precip': lambda t: precip_at_time_t,
                                 ...
                             }
                         }
                         If None, uses internal synthetic drivers.
    """
    start_time = time.time()
    
    print(f"=== Starting JULES Master Coupling: {weeks} weeks ===")
    
    # Debug: Check external drivers
    if external_drivers:
        print(f"External drivers provided: {list(external_drivers.keys())}")
        if 'era5_forcing' in external_drivers:
            print(f"ERA5 forcing variables: {list(external_drivers['era5_forcing'].keys())}")
        else:
            print("⚠️  WARNING: external_drivers provided but no 'era5_forcing' found!")
    else:
        print("Using internal synthetic forcing")
    
    # Initialize model states
    triffid_state = triffid_init if triffid_init is not None else [6.0, 3.0, 0.8, 0.2]
    rothc_state = rothc_init if rothc_init is not None else [C0_default[p] for p in POOLS]
    
    # Build default soil state if none provided
    if soil_init is None:
        nlay = SOIL_LAYERS['n_layers']
        # Create proper dict with both theta and T_soil as arrays
        soil_init = {
            'theta': np.full(nlay, INITIAL_CONDITIONS.get('theta_init', 0.3)),
            'T_soil': np.full(nlay, INITIAL_CONDITIONS.get('T_init', 283.15))
        }
    
    # Pre-allocate storage arrays for results
    results = {
        'triffid': np.zeros((4, weeks + 1)),
        'rothc': np.zeros((4, weeks + 1)),
        'weekly': {
            'nu_total': np.zeros(weeks + 1),
            'LAI_total': np.zeros(weeks + 1),
            'litterfall': np.zeros(weeks + 1),
            'soil_C_total': np.zeros(weeks + 1),
            'NPP': np.zeros(weeks),
            
            'Lb_tree': np.zeros(weeks + 1),
            'Lb_grass': np.zeros(weeks + 1), 
            'nu_tree': np.zeros(weeks + 1),
            'nu_grass': np.zeros(weeks + 1),
            
            'C_dpm': np.zeros(weeks + 1),
            'C_rpm': np.zeros(weeks + 1),
            'C_bio': np.zeros(weeks + 1),
            'C_hum': np.zeros(weeks + 1),
            
            'Pi_G': np.zeros(weeks),      
            'R_p': np.zeros(weeks),      
            'soil_moisture_mean': np.zeros(weeks),
            'soil_temp_mean': np.zeros(weeks),
        }
    }
    
    # Store initial conditions
    results['triffid'][:, 0] = triffid_state
    results['rothc'][:, 0] = rothc_state
    veg_vars = extract_vegetation_vars(triffid_state)
    
    # Initial detailed states
    results['weekly']['nu_total'][0] = veg_vars['nu_total']
    results['weekly']['LAI_total'][0] = veg_vars['LAI_total']
    results['weekly']['soil_C_total'][0] = sum(rothc_state)
    results['weekly']['Lb_tree'][0] = triffid_state[0]
    results['weekly']['Lb_grass'][0] = triffid_state[1]
    results['weekly']['nu_tree'][0] = triffid_state[2]
    results['weekly']['nu_grass'][0] = triffid_state[3]
    results['weekly']['C_dpm'][0] = rothc_state[0]
    results['weekly']['C_rpm'][0] = rothc_state[1]
    results['weekly']['C_bio'][0] = rothc_state[2]
    results['weekly']['C_hum'][0] = rothc_state[3]
    
    # Add archive for fast model results
    results['fast_results_archive'] = []
    
    # Store detailed weekly results
    results['weekly_results'] = []
    
    # Main simulation loop
    for week in range(weeks):
        # 1. Extract vegetation coupling variables from TRIFFID
        veg_vars = extract_vegetation_vars(triffid_state)
        
        # 2. Run fast models (EBM + Soil + Physiology)
        fast_results = run_fast_models_week(
            nu_cover=veg_vars['nu_total'],
            LAI_total=veg_vars['LAI_total'],
            soil_initial=soil_init,
            week_num=week,
            external_drivers=external_drivers
        )
        
        # Store the full fast_results (contains time series)
        results['weekly_results'].append(fast_results)
        
        # Archive first few weeks of detailed results for plotting
        if week < 4:  # Keep first 4 weeks of high-res data
            results['fast_results_archive'].append(fast_results.copy())
        
        # 3. Run slow models (TRIFFID + RothC)
        slow_results = run_slow_models_week(
            fast_results=fast_results,
            triffid_state=triffid_state,
            rothc_state=rothc_state,
            week_num=week
        )
        
        # 4. Update states for next week
        triffid_state = slow_results['triffid_new']
        rothc_state = slow_results['rothc_new']
        soil_state = fast_results['soil_final']
        
        # 5. Store weekly results
        results['triffid'][:, week+1] = triffid_state
        results['rothc'][:, week+1] = rothc_state
        
        results['weekly']['nu_total'][week+1] = slow_results['coupling_vars']['nu_total']
        results['weekly']['LAI_total'][week+1] = veg_vars['LAI_total']
        results['weekly']['litterfall'][week+1] = slow_results['coupling_vars']['litterfall']
        results['weekly']['soil_C_total'][week+1] = sum(rothc_state)
        results['weekly']['NPP'][week] = fast_results['NPP_total']
        
        results['weekly']['Lb_tree'][week+1] = triffid_state[0]
        results['weekly']['Lb_grass'][week+1] = triffid_state[1]
        results['weekly']['nu_tree'][week+1] = triffid_state[2]
        results['weekly']['nu_grass'][week+1] = triffid_state[3]
        
        results['weekly']['C_dpm'][week+1] = rothc_state[0]
        results['weekly']['C_rpm'][week+1] = rothc_state[1]
        results['weekly']['C_bio'][week+1] = rothc_state[2]
        results['weekly']['C_hum'][week+1] = rothc_state[3]
        
        results['weekly']['Pi_G'][week] = fast_results.get('Pi_G', 0)
        results['weekly']['R_p'][week] = fast_results.get('R_p', 0)
        results['weekly']['soil_moisture_mean'][week] = np.mean(fast_results['soil_theta'][0, :])
        results['weekly']['soil_temp_mean'][week] = np.mean(fast_results['soil_T'][0, :])
        
        print(f"Week {week+1:2d}: NPP={fast_results['NPP_total']:.4f} kg C/m²/wk, "
              f"θ={results['weekly']['soil_moisture_mean'][week]:.3f} m³/m³, "
              f"T={results['weekly']['soil_temp_mean'][week]-273.15:.1f}°C, "
              f"LAI={veg_vars['LAI_total']:.2f}")
    
    elapsed = time.time() - start_time
    print(f"\nSimulation completed in {elapsed:.1f} seconds")
    
    # Store metadata
    results['metadata'] = {
        'external_forcing': 'ERA5' if external_drivers else 'Synthetic',
        'weeks': weeks,
        'elapsed_time': elapsed
    }
    
    return results

# Legacy function for backward compatibility
def run_jules_simulation(config=None):
    """Legacy function - use jules_master_coupler directly"""
    if config is None:
        config = {'weeks': 52}
        
    results = jules_master_coupler(
        weeks=config.get('weeks', 52),
        triffid_init=config.get('triffid_init'),
        rothc_init=config.get('rothc_init'),
        soil_init=config.get('soil_init')
    )
    
    return results

if __name__ == '__main__':
    # Test with synthetic drivers
    results = jules_master_coupler(weeks=4)
    print(f"Test simulation complete. Used {results['metadata']['external_forcing']} forcing.")