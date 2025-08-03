# soil/simulation.py
import numpy as np
from scipy.integrate import solve_ivp

from julesf.soil.equations_moisture import moisture_rhs
from julesf.soil.equations_thermal import thermal_rhs
from julesf.soil.parameters import SOIL_LAYERS, SOIL_PROPERTIES, VAN_GENUCHTEN, THERMAL_PROPERTIES

def generate_soil_forcings(days=7, dt_hours=0.5):
    """
    Generate realistic soil forcing time series
    
    Parameters:
    - days:     simulation length (days)
    - dt_hours: time step (hours)
    
    Returns:
    - t:        time array (hours)
    - drivers:  dictionary of forcing functions
    """
    t = np.arange(0, days * 24, dt_hours)
    
    # Diurnal air temperature (K) - 10°C ± 5°C with random noise
    T_air = 283.15 + 5 * np.sin(2*np.pi*t/24) + np.random.normal(0, 1, len(t))
    
    # Precipitation (kg/m²/s) - random bursts, positive half-normal
    precip_base = np.maximum(0, np.random.normal(0, 1e-5, len(t)))
    precip_events = (np.random.random(len(t)) < 0.02)  # 2% chance of rain each timestep
    precip = precip_base + precip_events * np.random.exponential(5e-5, len(t))
    
    # Evapotranspiration (kg/m²/s) - daytime peak, temperature dependent
    ET_potential = np.maximum(0, 2e-5 * np.sin(np.pi * (t % 24) / 24))
    ET = ET_potential * np.maximum(0, (T_air - 273.15) / 20)  # Temperature scaling
    
    # Surface heat flux (W/m²) - diurnal solar cycle
    G_surface = 100 * np.sin(2*np.pi*t/24) * np.maximum(0, np.sin(np.pi * (t % 24) / 24))
    
    # Create interpolation functions for drivers
    from scipy.interpolate import interp1d
    drivers = {
        'air_temperature': interp1d(t, T_air, bounds_error=False, fill_value='extrapolate'),
        'precipitation': interp1d(t, precip, bounds_error=False, fill_value=0),
        'evapotranspiration': interp1d(t, ET, bounds_error=False, fill_value=0),
        'surface_heat_flux': interp1d(t, G_surface, bounds_error=False, fill_value=0),
    }
    
    return t, drivers

def soil_rhs(t, y, params, drivers):
    """
    Combined RHS for coupled soil moisture-thermal system
    
    State vector: y = [θ₁, θ₂, θ₃, θ₄, T₁, T₂, T₃, T₄]
    
    Parameters:
    - t:       current time (hours)
    - y:       state vector
    - params:  soil parameter dictionary  
    - drivers: forcing functions
    
    Returns:
    - dydt: state derivative vector
    """
    n_layers = params['n_layers']
    
    theta = y[:n_layers]           # moisture content (m³/m³)
    T_soil = y[n_layers:]          # soil temperature (K)
    
    t_seconds = t * 3600
    dtheta_dt = moisture_rhs(theta, T_soil, params, drivers, t)
    
    from .equations_moisture import vertical_water_flux, infiltration_rate
    W_flux = vertical_water_flux(theta, params['layer_thickness'], params)
    W_flux[0] = infiltration_rate(params, drivers, t)
    
    W_layer = np.zeros(n_layers)
    for k in range(n_layers):
        W_layer[k] = 0.5 * (W_flux[k] + W_flux[k+1])
    
    dT_dt = thermal_rhs(T_soil, theta, W_layer, params, drivers, t)
    
    dtheta_dt = dtheta_dt * 3600  # m³/m³/hr
    dT_dt = dT_dt * 3600          # K/hr
    
    dydt = np.concatenate([dtheta_dt, dT_dt])
    
    return dydt

def solve_soil_rhs(t_span, theta_init, T_init, drivers, method="BDF", dt_out=0.5):
    """
    IVP solver for coupled soil moisture-thermal system
    
    Parameters:
    - t_span:     time span (start, end) in hours
    - theta_init: initial moisture content array (m³/m³)
    - T_init:     initial temperature array (K)  
    - drivers:    forcing functions dictionary
    - method:     ODE solver method
    - dt_out:     output timestep (hours)
    
    Returns:
    - t:       time array (hours)
    - theta:   moisture content array (n_layers × n_times)
    - T_soil:  soil temperature array (n_layers × n_times)
    """
    
    # Combine soil parameters
    params = {**SOIL_LAYERS, **SOIL_PROPERTIES, **VAN_GENUCHTEN, **THERMAL_PROPERTIES}
    
    # Set up initial state vector
    y0 = np.concatenate([theta_init, T_init])
    
    # Time evaluation points
    t_eval = np.arange(t_span[0], t_span[1] + dt_out, dt_out)
    
    # Solve ODE system 
    print(f"Solving coupled soil system from t={t_span[0]} to t={t_span[1]} hours...")
    print(f"State variables: {len(y0)} ({params['n_layers']} moisture + {params['n_layers']} temperature)")
    
    sol = solve_ivp(
        fun=lambda t, y: soil_rhs(t, y, params, drivers),
        t_span=t_span,
        y0=y0,
        method=method,
        t_eval=t_eval,
        rtol=1e-6,
        atol=1e-9
    )
    
    if not sol.success:
        raise RuntimeError(f"ODE solver failed: {sol.message}")
    
    # Split solution
    n_layers = params['n_layers']
    theta = sol.y[:n_layers, :]      # moisture (n_layers × n_times)
    T_soil = sol.y[n_layers:, :]     # temperature (n_layers × n_times)
    
    print(f"Simulation completed successfully!")
    print(f"Final moisture range: {theta.min():.3f} - {theta.max():.3f} m³/m³")
    print(f"Final temperature range: {T_soil.min():.1f} - {T_soil.max():.1f} K")
    
    return sol.t, theta, T_soil
