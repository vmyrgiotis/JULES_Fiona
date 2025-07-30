# TRIFFID parameters (Clark et al 2011, JULES standard units)
import numpy as np

PFT_PARAMS = {
    'broadleaf_tree': {
        'sigma1': 0.0375,       # kg C m⁻² LAI⁻¹, leaf carbon per LAI
        'awl': 0.65,            # kg C m⁻², woody allometry coefficient
        'bwl': 1.667,           # dimensionless, woody allometry exponent
        'gamma_l': 0.25,        # yr⁻¹, leaf turnover rate
        'gamma_r': 0.25,        # yr⁻¹, root turnover rate  
        'gamma_w': 0.005,       # yr⁻¹, wood turnover rate
        'Lmin': 1.0,            # m² m⁻², minimum LAI
        'Lmax': 9.0,            # m² m⁻², maximum LAI
    },
    'C3_grass': {
        'sigma1': 0.0250,
        'awl': 0.005,
        'bwl': 1.667,
        'gamma_l': 0.25,
        'gamma_r': 0.25,
        'gamma_w': 0.20,
        'Lmin': 1.0,
        'Lmax': 4.0,
    }
}

COMPETITION_MATRIX = {
    'c': np.array([[1.0, 1.0], [0.0, 1.0]]),  # Competition coefficients
}

INITIAL_CONDITIONS = {
    'Lb_init': np.array([6.0, 3.0]),    # m² m⁻², initial balanced LAI
    'nu_init': np.array([0.8, 0.2]),    # dimensionless, initial fractional cover
}