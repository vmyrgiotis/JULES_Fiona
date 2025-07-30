# soil/parameters.py

import numpy as np

# Layer configuration (from JULES paper)
SOIL_LAYERS = {
    'n_layers': 4,
    'depths': np.array([0.1, 0.25, 0.65, 2.0]),      # m, cumulative depth
    'layer_thickness': np.array([0.1, 0.15, 0.4, 1.35]),   # m, ΔZₖ layer thickness - FIXED KEY NAME
}

# Soil physical properties (LOAM soil type) 
# TODO: Research exact JULES default values
SOIL_PROPERTIES = {
    'theta_sat': 0.43,      # m³/m³, saturated moisture content (θₛ)
    'theta_res': 0.078,     # m³/m³, residual moisture content (θᵣ) 
    'K_sat': 1.16e-6,       # m/s, saturated hydraulic conductivity (Kₕₛ)
    'rho_soil': 1300.0,     # kg/m³, soil bulk density
    'c_soil': 800.0,        # J/kg/K, soil specific heat capacity
}

# van Genuchten parameters (Equations 60-61)
# TODO: Find JULES-specific values for different soil types
VAN_GENUCHTEN = {
    'alpha_vG': 3.6,        # m⁻¹, inverse air entry potential
    'n_vG': 1.56,           # dimensionless, pore size distribution
    'm_vG': 0.359,          # = 1 - 1/n_vG
    'l': 0.5,               # dimensionless, tortuosity parameter
}

# Thermal properties (Equations 62-67)
# TODO: Verify these values against JULES documentation
THERMAL_PROPERTIES = {
    'lambda_dry': 0.25,     # W/m/K, dry soil thermal conductivity (λdry)
    'lambda_sat': 2.2,      # W/m/K, saturated soil thermal conductivity (λₐ)
    'c_water': 4180.0,      # J/kg/K, specific heat of water (Cwater)
    'rho_water': 1000.0,    # kg/m³, density of water
}

# Initial conditions
INITIAL_CONDITIONS = {
    'theta_init': 0.3,      # m³/m³, field capacity for all layers
    'T_init': 283.15,       # K, initial soil temperature
}

# Time stepping
TIME_SETTINGS = {
    'dt_hours': 0.5,        # hours, timestep to match EBM/Physiology
    'days': 7,              # days, simulation length
}