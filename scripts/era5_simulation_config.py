import os

# Data file paths
DATA_DIR = "../data"
ERA5_DATA_FILE = os.path.join(DATA_DIR, "era5_jules_met_ins_2002-2016.csv")
ERA5_METADATA_FILE = os.path.join(DATA_DIR, "era5_jules_met_ins_metadata.csv")

# Simulation configurations
SIMULATION_CONFIGS = {
    'short_test': {
        'weeks': 4,
        'description': 'Short 4-week test simulation',
        'output_dir': 'results/era5_test'
    },
    
    'seasonal': {
        'weeks': 13,
        'description': 'Seasonal simulation (1 quarter)',
        'output_dir': 'results/era5_seasonal'
    },
    
    'annual': {
        'weeks': 52,
        'description': 'Full annual simulation',
        'output_dir': 'results/era5_annual'
    },
    
    'multi_year': {
        'weeks': 15 * 52,      # ~15 years (2002–2016)
        'description': '15-year simulation (2002–2016)',
        'output_dir': 'results/era5_multi_year'
    }
}

# Model initial condition
INITIAL_CONDITIONS = {
    'default': {
        'triffid_init': None,  # Use model defaults
        'rothc_init': None,    # Use model defaults
        'soil_init': None      # Use model defaults
    },
    
    'custom': {
        'triffid_init': [5.0, 4.0, 0.7, 0.3],  # [Lb_tree, Lb_grass, nu_tree, nu_grass]
        'rothc_init': None,    # Could specify carbon pools
        'soil_init': None      # Could specify soil state
    }
}

def get_config(config_name='annual', initial_conditions='default'):
    if config_name not in SIMULATION_CONFIGS:
        raise ValueError(f"Unknown config: {config_name}. Available: {list(SIMULATION_CONFIGS.keys())}")
    
    if initial_conditions not in INITIAL_CONDITIONS:
        raise ValueError(f"Unknown initial conditions: {initial_conditions}. Available: {list(INITIAL_CONDITIONS.keys())}")
    
    config = SIMULATION_CONFIGS[config_name].copy()
    config.update(INITIAL_CONDITIONS[initial_conditions])
    config['data_file'] = ERA5_DATA_FILE
    config['metadata_file'] = ERA5_METADATA_FILE
    
    return config