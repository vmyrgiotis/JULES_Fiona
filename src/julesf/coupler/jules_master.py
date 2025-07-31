import numpy as np
import time

from julesf.coupler.fast_models_wrapper import run_fast_models_week
from julesf.coupler.slow_models_wrapper import run_slow_models_week
from julesf.coupler.coupling_utils import extract_vegetation_vars
from julesf.rothc.parameters import POOLS, C0_default
from julesf.soil.parameters import SOIL_LAYERS, INITIAL_CONDITIONS

def jules_master_coupler(weeks=52, triffid_init=None, rothc_init=None, soil_init=None):
    """
    Run coupled JULES model for specified number of weeks
    
    Args:
        weeks: Number of weeks to run
        triffid_init: Optional initial TRIFFID state [Lb_tree, Lb_grass, nu_tree, nu_grass]
        rothc_init: Optional initial RothC state [C_dpm, C_rpm, C_bio, C_hum]
        soil_init: Optional initial soil state {'theta': [...], 'T_soil': [...]}
    
    Returns:
        Dict with full simulation results
    """
    start_time = time.time()
    
    # Initialize model states
    triffid_state = triffid_init if triffid_init is not None else [6.0, 3.0, 0.8, 0.2]
    rothc_state = rothc_init if rothc_init is not None else [C0_default[p] for p in POOLS]
    
    # Build default soil state if none provided
    if soil_init is None:
        n = SOIL_LAYERS['n_layers']
        soil_state = {
            'theta': np.full(n, INITIAL_CONDITIONS['theta_init']),
            'T_soil': np.full(n, INITIAL_CONDITIONS['T_init'])
        }
    else:
        soil_state = soil_init
    
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
        }
    }
    
    # Store initial conditions
    results['triffid'][:, 0] = triffid_state
    results['rothc'][:, 0] = rothc_state
    veg_vars = extract_vegetation_vars(triffid_state)
    results['weekly']['nu_total'][0] = veg_vars['nu_total']
    results['weekly']['LAI_total'][0] = veg_vars['LAI_total']
    results['weekly']['soil_C_total'][0] = sum(rothc_state)
    
    print(f"=== Starting JULES Master Coupling: {weeks} weeks ===")
    
    for week in range(weeks):
        print(f"Week {week+1}/{weeks}", end="\r")
        
        # 1. Extract vegetation coupling variables
        veg_vars = extract_vegetation_vars(triffid_state)
        
        # 2. Run fast models (0.5h) for 1 week
        fast_results = run_fast_models_week(
            nu_cover=veg_vars['nu_total'],
            LAI_total=veg_vars['LAI_total'],
            soil_initial=soil_state,
            week_num=week
        )
        
        # 3. Run slow models (weekly)
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
    
    elapsed = time.time() - start_time
    print(f"\nSimulation completed in {elapsed:.1f} seconds")
    
    return results

def run_jules_simulation(config=None):
    """User-friendly interface with configuration options"""
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
    results = run_jules_simulation()