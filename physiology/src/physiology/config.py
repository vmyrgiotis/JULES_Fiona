# PFT parameters
PFT_PARAMS = {
    'C3_grass': {
        'alpha': 0.12,
        'omega': 0.15,
        'f_dr': 0.015,
        'n0': 0.073,
        'n_e': 0.0008,
        'Q10_Kc': 2.1,
        'Q10_Ko': 1.2,
        'Q10_rs': 0.57,
        'T_low': 0.0,
        'T_upp': 36.0,
    },
    # add other PFTs here
}

# Simulation settings
SIM_SETTINGS = {
    'days': 10,
    'dt_hours': 1,
    'T_mean': 20.0,
    'T_amp': 5.0,
    'I_max': 1500.0,
    'ca_ppm': 400.0,
    'chi': 0.7,
    'P': 101325.0,
    'O2_fraction': 0.21,
}
