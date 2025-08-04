# src/soil/equations_moisture.py
"""
JULES Soil Moisture Equations
Main equation: Eq. (56) - Richards equation for soil water flow
"""
import numpy as np

from julesf.soil.hydraulic_options import simple_hydraulic_conductivity, simple_matric_potential

def vertical_water_flux(theta, layer_thickness, params):
    """
    Vertical water flux between layers - Wₖ
    JULES Equation (57): W = Kₕ(∂ψ/∂z + 1)
    
    Parameters:
    - theta:           moisture content array (m³/m³)  
    - layer_thickness: ΔZₖ array (m)
    - params:          soil parameters dict
    
    Returns:
    - W: water flux array (kg/m²/s), positive = downward
    """
    n_layers = len(theta)
    W = np.zeros(n_layers + 1)  # fluxes at layer interfaces (0 to n_layers)
    
    Kh = np.array([simple_hydraulic_conductivity(theta[k], params) 
                   for k in range(n_layers)])
    psi = np.array([simple_matric_potential(theta[k], params) 
                    for k in range(n_layers)])
    
    W[0] = 0.0  # Surface flux W₀ (boundary condition), will be set by infiltration from drivers
    
    for k in range(1, n_layers):
        # Average hydraulic conductivity between layers
        Kh_avg = 0.5 * (Kh[k-1] + Kh[k])
        
        # Matric potential gradient ∂ψ/∂z (finite difference)
        dz = 0.5 * (layer_thickness[k-1] + layer_thickness[k])
        dpsi_dz = (psi[k] - psi[k-1]) / dz
        
        # Darcy's law: W = Kₕ(∂ψ/∂z + 1) 
        W[k] = Kh_avg * (dpsi_dz + 1.0)
    
    W[n_layers] = Kh[-1]  # Bottom boundary, gravitational drainage from bottom layer
    
    return W

def evapotranspiration_extraction(theta, params, drivers, t):
    """
    Evapotranspiration extraction from each layer - Eₖ  
    
    Parameters:
    - theta:   moisture content array (m³/m³)
    - params:  soil parameters dict
    - drivers: forcing functions
    - t:       current time
    
    Returns:
    - E: evapotranspiration array (kg/m²/s)
    """
    n_layers = len(theta)
    E = np.zeros(n_layers)
    
    if 'evapotranspiration' in drivers:
        ET_total = drivers['evapotranspiration'](t)  # kg/m²/s
        
        # Distribute ET exponentially with depth (most from surface layers)
        layer_weights = np.exp(-np.arange(n_layers))  # exponential decay
        layer_weights = layer_weights / np.sum(layer_weights)  # normalize
        
        for k in range(n_layers):
            # Extract proportional to layer weight and available water
            theta_available = max(0, theta[k] - params['theta_res'])
            if theta_available > 0:
                E[k] = ET_total * layer_weights[k]
            else:
                E[k] = 0.0  # No extraction if below residual
    
    return E

def lateral_runoff(theta, params):
    """
    Lateral runoff from each layer - Rbk
    Set to zero for standalone soil column (no TOPMODEL)
    
    Parameters:
    - theta:  moisture content array (m³/m³)
    - params: soil parameters dict
    
    Returns:
    - R: lateral runoff array (kg/m²/s) - all zeros for now
    """
    n_layers = len(theta)
    R = np.zeros(n_layers)
    
    # No lateral runoff for standalone soil column
    # TODO: Implement TOPMODEL Rbk for catchment coupling
    
    return R

def infiltration_rate(params, drivers, t):
    """
    Surface infiltration rate - W₀
    JULES Equation (49): W₀ = Σⱼ νⱼ(Tⱼ + Sₘⱼ - Yⱼ)
    
    Simplified: W₀ = max(0, precipitation - surface_runoff)
    
    Parameters:
    - params:  soil parameters dict
    - drivers: forcing functions  
    - t:       current time
      
    Returns:
    - W0: infiltration rate (kg/m²/s)
    """
    if 'precipitation' in drivers:
        precip = drivers['precipitation'](t)  # kg/m²/s
        
        # Convert K_sat from m/s to kg/m²/s properly
        rho_water = params.get('rho_water', 1000.0)  # kg/m³
        K_sat_kg = params['K_sat'] * rho_water  # m/s × kg/m³ = kg/m²/s
        
        W0 = min(precip, K_sat_kg)
        return max(0, W0)
    else:
        return 0.0

def moisture_rhs(theta, T_soil, params, drivers, t):
    """
    Right-hand side of moisture equation - JULES Equation (56)
    
    dθₖ/dt = Wₖ₋₁ - Wₖ - Eₖ - Rbk
    
    For standalone soil: Rbk = 0, so dθₖ/dt = Wₖ₋₁ - Wₖ - Eₖ
    """
    layer_thickness = params['layer_thickness']
    n_layers = len(theta)
    
    W = vertical_water_flux(theta, layer_thickness, params)    
    W[0] = infiltration_rate(params, drivers, t)
    
    E = evapotranspiration_extraction(theta, params, drivers, t)
    R = lateral_runoff(theta, params) 
    
    dtheta_dt = np.zeros(n_layers)
    
    for k in range(n_layers):
        water_flux_conv = W[k] - W[k+1] - E[k] - R[k]  
        rho_water = params['rho_water']  # Convert from kg/m²/s to m³/m³/s
        dtheta_dt[k] = water_flux_conv / (layer_thickness[k] * rho_water)
    
    return dtheta_dt
