#parameters.py
PHOTOSYNTHESIS_PFT_PARAMS = {
    'broadleaf_tree': {
        'alpha':   0.08,
        'omega':   0.15,
        'f_dr':    0.015,
        'n0':      0.046,
        'n_e':     0.0008,
        'Q10_Kc':  2.1,
        'Q10_Ko':  1.2,
        'Q10_rs':  0.57,
        'T_low':   0.0,
        'T_upp':   36.0,
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
        'T_low':  -10.0,
        'T_upp':   26.0,
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
        'T_low':   0.0,
        'T_upp':   36.0,
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
        'T_low':  13.0,
        'T_upp':   45.0,
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
        'T_low':   0.0,
        'T_upp':   36.0,
    },
}


CANOPY_PARAMS = {
    'k_ext':        0.5,
    'LAI':          3.0,
    'h':            1
}

RESPIRATION_PARAMS = {
    'beta':         1.0,
    'rg':           0.25,         # growth resp coeff (25%)
    'n_m':           0.2,

    'sigma1':    0.015,           # σ1₁
    'mean_leaf_N':        0.033,  # nₘ: mean leaf nitrogen concentration, kg N per kg C
    'mu_root_leaf_N':     1.0,    # μ₍rl₎: root N / leaf N
    'mu_stem_leaf_N':     0.1,    # μ₍sl₎: stem N / leaf N
    'eta_root_C': 0.20,           # η₍rl_C₎: kg C root per unit LAI
    'eta_stem_C':         0.02,   # η₍sl_C₎: kg C stem per m height per unit LAI
}

SIM_SETTINGS = {
    'days':         10,
    'dt_hours':     1,
    'T_mean':       20.0,
    'T_amp':        5.0,
    'I_max':        1500.0,
    'ca_ppm':       400.0,
    'chi':          0.7,
    'P':            101325.0,
    'O2_fraction':  0.21
}