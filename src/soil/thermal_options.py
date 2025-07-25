# src/soil/thermal_options.py
import numpy as np

def simple_thermal_conductivity(theta, params):
    """
    Here, I am using a simplified approach to calculate thermal conductivity.
    TODO: Implement Equations (62-67) with Kersten number
    
    Parameters:
    - theta: volumetric moisture content (m³/m³)
    
    Returns:
    - lambda_soil (λ): thermal conductivity (W/m/K)
    """
    lambda_dry = params['lambda_dry'] 
    lambda_sat = params['lambda_sat']
    theta_sat = params['theta_sat']
    
    Se = np.clip(theta / theta_sat, 0, 1)
    lambda_soil = lambda_dry + (lambda_sat - lambda_dry) * Se
    
    return lambda_soil