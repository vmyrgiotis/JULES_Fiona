# src/soil/hydraulic_options.py
import numpy as np

def simple_hydraulic_conductivity(theta, params):
    """
    Here, I am using a simplified approach to calculate hydraulic conductivity.
    TODO: Implement full Equation (61): Kₕ = Kₕₛ Sᵉˡ [1-(1-Sᵉ^(1/m))^m]²
    
    Parameters:
    - theta: volumetric moisture content (m³/m³)
    - params: soil parameters dict
    
    Returns:
    - Kh: hydraulic conductivity (m/s)
    """
    theta_sat = params['theta_sat']
    K_sat = params['K_sat']
    
    Se = np.clip(theta / theta_sat, 0, 1)  
    Kh = K_sat * Se**3
    
    return Kh

def simple_matric_potential(theta, params):
    """
    Now, I am using a simplified approach to calculate matric potential
    TODO: Implement full Equation (60): θ = θᵣ + (θₛ-θᵣ)/[1+(αψ)ⁿ]^m
    
    Parameters:
    - theta: volumetric moisture content (m³/m³)
    
    Returns:
    - psi: matric potential (m, negative values)
    """
    theta_sat = params['theta_sat']
    theta_res = params['theta_res']

    Se = (theta - theta_res) / (theta_sat - theta_res)
    Se = np.clip(Se, 0.01, 1.0)
    psi = -1.0 / Se**2 
    
    return psi