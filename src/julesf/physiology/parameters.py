# parameters.py 

PHOTOSYNTHESIS_PFT_PARAMS = {
    'broadleaf_tree': {
        'alpha':   0.08,          # mol CO₂ mol⁻¹ photon, quantum efficiency
        'omega':   0.15,          # dimensionless, leaf scattering coefficient
        'f_dr':    0.015,         # dimensionless, dark respiration coefficient
        'n0':      0.046,         # kg N kg⁻¹ C, top-leaf N concentration
        'n_e':     0.0008,        # mol CO₂ m⁻² s⁻¹ per kg C, conversion factor
        'Q10_Kc':  2.1,           # dimensionless, Q10 for Kc
        'Q10_Ko':  1.2,           # dimensionless, Q10 for Ko
        'Q10_rs':  0.57,          # dimensionless, Q10 for RuBisCO specificity
        'T_low':   273.15,        # K, lower temperature parameter (0°C)
        'T_upp':   309.15,        # K, upper temperature parameter (36°C)
    },
    'needleleaf_tree': {
        'alpha':   0.08,
        'omega':   0.15, 
        'f_dr':    0.015,
        'n0':      0.033,
        'n_e':     0.0008,
        'Q10_Kc':  2.1,
        'Q10_Ko':  1.2,
        'Q10_rs':  0.57,
        'T_low':   263.15,        # K, (-10°C)
        'T_upp':   299.15,        # K, (26°C)
    },
    'C3_grass': {
        'alpha':   0.12,
        'omega':   0.15,
        'f_dr':    0.015,
        'n0':      0.073,
        'n_e':     0.0008,
        'Q10_Kc':  2.1,
        'Q10_Ko':  1.2,
        'Q10_rs':  0.57,
        'T_low':   273.15,        # K, (0°C)
        'T_upp':   309.15,        # K, (36°C)
    },
    'C4_grass': {
        'alpha':   0.06,
        'omega':   0.17,
        'f_dr':    0.025,
        'n0':      0.060,
        'n_e':     0.0004,
        'Q10_Kc':  2.1,
        'Q10_Ko':  1.2,
        'Q10_rs':  0.57,
        'T_low':   286.15,        # K, (13°C)
        'T_upp':   318.15,        # K, (45°C)
    },
    'shrub': {
        'alpha':   0.08,
        'omega':   0.15,
        'f_dr':    0.015,
        'n0':      0.060,
        'n_e':     0.0008,
        'Q10_Kc':  2.1,
        'Q10_Ko':  1.2,
        'Q10_rs':  0.57,
        'T_low':   273.15,        # K, (0°C)
        'T_upp':   309.15,        # K, (36°C)
    },
}

CANOPY_PARAMS = {
    'k_ext':        0.5,          # dimensionless, light extinction coefficient
    'LAI':          3.0,          # m² m⁻², leaf area index
    'h':            1.0           # m, vegetation height
}

RESPIRATION_PARAMS = {
    'beta':         1.0,          # dimensionless, temperature response parameter
    'rg':           0.25,         # dimensionless, growth respiration coefficient (25%)
    'n_m':          0.2,          # dimensionless, maintenance respiration parameter
    
    'sigma1':       0.015,        # kg C m⁻² (LAI)⁻¹, specific leaf density (σ₁₁)
    'mean_leaf_N':  0.033,        # kg N kg⁻¹ C, mean leaf nitrogen concentration (nₘ)
    'mu_root_leaf_N': 1.0,        # dimensionless, root N / leaf N ratio (μ₍rl₎)
    'mu_stem_leaf_N': 0.1,        # dimensionless, stem N / leaf N ratio (μ₍sl₎)
    'eta_root_C':   0.20,         # kg C root per unit LAI (η₍rl_C₎)
    'eta_stem_C':   0.02,         # kg C stem m⁻¹ height per unit LAI (η₍sl_C₎)
}

SIM_SETTINGS = {
    'days':         10,
    'dt_hours':     1,
    'T_mean':       293.15,       # K, mean temperature (20°C)
    'T_amp':        5.0,          # K, temperature amplitude 
    'I_max':        1500.0,       # μmol m⁻² s⁻¹, maximum PAR
    'ca_ppm':       400.0,        # ppm, ambient CO₂ concentration  
    'chi':          0.7,          # dimensionless, ci/ca ratio
    'P':            101325.0,     # Pa, atmospheric pressure
    'O2_fraction':  0.21          # dimensionless, O₂ fraction
}