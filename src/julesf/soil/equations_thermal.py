# src/soil/equations_thermal.py
"""
JULES Soil Thermal Equations
Main equation: Eq. (53) - Heat conduction with advection
"""
import numpy as np

from julesf.soil.thermal_options import simple_thermal_conductivity

def volumetric_heat_capacity(theta, params):
    """
    Volumetric heat capacity of soil-water mixture - Cₖ
    Used in Equation (53): Cₖ ΔZₖ dTₖ/dt = Gₖ₋₁ - Gₖ - Jₖ ΔZₖ
    
    Parameters:
    - theta:   volumetric moisture content (m³/m³)
    - params:  soil parameters dict
    
    Returns:
    - C:       volumetric heat capacity (J/m³/K)
    """
    rho_soil = params['rho_soil']
    c_soil = params['c_soil'] 
    rho_water = params['rho_water']
    c_water = params['c_water']
    
    C = rho_soil * c_soil + theta * rho_water * c_water
    
    return C

def diffusive_heat_flux(T_soil, layer_thickness, params, theta):
    """
    Diffusive heat flux between layers - Gₖ
    JULES Equation (54): G = λₛₒᵢₗ ∂T/∂z
    
    Parameters:
    - T_soil:          temperature array for all layers (K)
    - layer_thickness: ΔZₖ array (m)
    - params:          soil parameters dict
    - theta:           moisture content array (m³/m³)
    
    Returns:
    - G: heat flux array (W/m²), positive = downward
    """
    n_layers = len(T_soil)
    G = np.zeros(n_layers + 1)  # fluxes at layer interfaces (0 to n_layers)
    
    lambda_soil = np.array([simple_thermal_conductivity(theta[k], params) 
                           for k in range(n_layers)])

    G[0] = 0.0  # Surface flux (boundary condition) - will be set by surface energy balance

    for k in range(1, n_layers):
        # Average thermal conductivity between layers
        lambda_avg = 0.5 * (lambda_soil[k-1] + lambda_soil[k])
        
        # Temperature gradient ∂T/∂z (finite difference)
        dz = 0.5 * (layer_thickness[k-1] + layer_thickness[k])
        dT_dz = (T_soil[k] - T_soil[k-1]) / dz
        
        # Fourier's law: G = -λ ∂T/∂z (negative sign for downward positive)
        G[k] = lambda_avg * dT_dz
    
    # Bottom boundary - free heat conduction
    lambda_bottom = lambda_soil[-1]
    # Assume zero gradient at bottom for now
    G[n_layers] = 0.0
    
    return G

def advective_heat_flux(T_soil, W_flux, layer_thickness, params):
    """
    Advective heat flux by water movement - Jₖ
    JULES Equation (55): J = C_water W ∂T/∂z
    
    Parameters:
    - T_soil:          temperature array (K)
    - W_flux:          water flux array from moisture equation (kg/m²/s)
    - layer_thickness: ΔZₖ array (m)
    - params:          soil parameters dict
    
    Returns:
    - J: advective heat flux (W/m²)
    """
    n_layers = len(T_soil)
    J = np.zeros(n_layers)
    
    c_water = params['c_water']  # J/kg/K
    
    for k in range(n_layers):
        if k == 0:
            # Surface layer - use surface temperature gradient
            if n_layers > 1:
                dz = 0.5 * (layer_thickness[0] + layer_thickness[1])
                dT_dz = (T_soil[1] - T_soil[0]) / dz
            else:
                dT_dz = 0.0
        elif k == n_layers - 1:
            # Bottom layer - use bottom temperature gradient
            dz = 0.5 * (layer_thickness[k-1] + layer_thickness[k])
            dT_dz = (T_soil[k] - T_soil[k-1]) / dz
        else:
            # Middle layers - central difference
            dz = 0.5 * (layer_thickness[k-1] + layer_thickness[k+1])
            dT_dz = (T_soil[k+1] - T_soil[k-1]) / dz
        
        # Advective flux: J = C_water * W * ∂T/∂z
        J[k] = c_water * W_flux[k] * dT_dz
    
    return J

def thermal_rhs(T_soil, theta, W_flux, params, drivers, t):
    """
    Right-hand side of thermal equation - JULES Equation (53)
    
    Cₖ ΔZₖ dTₖ/dt = Gₖ₋₁ - Gₖ - Jₖ ΔZₖ
    
    Rearranged: dTₖ/dt = (Gₖ₋₁ - Gₖ - Jₖ ΔZₖ) / (Cₖ ΔZₖ)
    
    Parameters:
    - T_soil:    soil temperature array (K)
    - theta:     moisture content array (m³/m³)
    - W_flux:    water flux array (kg/m²/s)
    - params:    soil parameters dict
    - drivers:   forcing functions
    - t:         current time
    
    Returns:
    - dT_dt: temperature change rate (K/s)
    """
    layer_thickness = params['layer_thickness']
    n_layers = len(T_soil)
    
    C = np.array([volumetric_heat_capacity(theta[k], params) 
                  for k in range(n_layers)])
    
    G = diffusive_heat_flux(T_soil, layer_thickness, params, theta)
    J = advective_heat_flux(T_soil, W_flux, layer_thickness, params)
    dT_dt = np.zeros(n_layers)
    
    for k in range(n_layers):
        heat_flux_conv = G[k] - G[k+1] - J[k] * layer_thickness[k]        
        dT_dt[k] = heat_flux_conv / (C[k] * layer_thickness[k])
    
    # Surface boundary condition 
    if 'surface_heat_flux' in drivers:
        G_surface = drivers['surface_heat_flux'](t)  # W/m²
        dT_dt[0] += G_surface / (C[0] * layer_thickness[0])
    
    return dT_dt
