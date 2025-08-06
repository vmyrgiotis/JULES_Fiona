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
    
    # Base precipitation: very light continuous
    precip_base = np.maximum(0, np.random.normal(0, 1e-4, len(t)))  
    
    # Rain events: 2% chance per timestep, with realistic intensity
    precip_events = (np.random.random(len(t)) < 0.02)  # 2% chance of rain
    rain_intensity = np.random.exponential(1e-3, len(t)) 
    precip = precip_base + precip_events * rain_intensity
    
    # Evapotranspiration (kg/m²/s) 
    # Typical ET: 2-6 mm/day = 2.3e-5 to 6.9e-5 kg/m²/s average
    ET_potential = np.maximum(0, 1e-4 * np.sin(np.pi * (t % 24) / 24))  
    ET = ET_potential * np.maximum(0, (T_air - 273.15) / 20)  # Temperature scaling
    
    # Surface heat flux (W/m²) - diurnal solar cycle
    G_surface = 100 * np.sin(2*np.pi*t/24) * np.maximum(0, np.sin(np.pi * (t % 24) / 24))
    
    # Create interpolation functions for drivers
    from scipy.interpolate import interp1d

    rho_water = 1000.0  # kg/m³

    # Convert precip depth (m per dt_hours) → flux kg/m²/s
    precip_flux = precip * rho_water / (dt_hours * 3600.0)

    drivers = {
        'air_temperature': interp1d(t, T_air, bounds_error=False, fill_value='extrapolate'),
        'precipitation':   interp1d(t, precip_flux, bounds_error=False, fill_value=0.0),
        'evapotranspiration': interp1d(t, ET, bounds_error=False, fill_value=0.0),
        'surface_heat_flux': interp1d(t, G_surface, bounds_error=False, fill_value=0.0),
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
    """
    
    # Combine soil parameters
    params = {**SOIL_LAYERS, **SOIL_PROPERTIES, **VAN_GENUCHTEN, **THERMAL_PROPERTIES}
    
    # FIXED: Ensure both theta_init and T_init are arrays with same length
    if theta_init is None:
        n_layers = params['n_layers']
        theta_init = np.full(n_layers, INITIAL_CONDITIONS.get('theta_init', 0.3))
    
    if T_init is None:
        n_layers = len(theta_init)  # Use theta_init length
        T_init = np.full(n_layers, INITIAL_CONDITIONS.get('T_init', 283.15))
    
    # Ensure T_init is an array, not a scalar
    if np.isscalar(T_init):
        n_layers = len(theta_init)
        T_init = np.full(n_layers, T_init)
    
    # Ensure both arrays have the same length
    if len(theta_init) != len(T_init):
        n_layers = params['n_layers']
        print(f"⚠️  Warning: theta_init ({len(theta_init)}) and T_init ({len(T_init)}) have different lengths")
        print(f"   Resizing both to {n_layers} layers")
        theta_init = np.full(n_layers, np.mean(theta_init) if len(theta_init) > 0 else 0.3)
        T_init = np.full(n_layers, np.mean(T_init) if len(T_init) > 0 else 283.15)
    
    # Set up initial state vector
    y0 = np.concatenate([theta_init, T_init])
    
    # Time evaluation points
    t_eval = np.arange(t_span[0], t_span[1] + dt_out, dt_out)
    
    # Solve ODE system 
    print(f"Solving coupled soil system from t={t_span[0]} to t={t_span[1]} hours...")
    print(f"State variables: {len(y0)} ({len(theta_init)} moisture + {len(T_init)} temperature)")
    
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
    n_layers = len(theta_init)
    theta    = sol.y[:n_layers, :]      # moisture (n_layers × n_times)
    T_soil   = sol.y[n_layers:, :]      # temperature (n_layers × n_times)
 
    # ENFORCE physical bounds on θ
    θ_res = params.get('theta_res', 0.0)
    θ_sat = params.get('theta_sat', 1.0)
    theta = np.clip(theta, θ_res, θ_sat)
 
    print(f"Simulation completed successfully!")
    print(f"Final moisture range (clamped): {theta.min():.3f} - {theta.max():.3f} m³/m³")
    print(f"Final temperature range: {T_soil.min():.1f} - {T_soil.max():.1f} K")
 
    return sol.t, theta, T_soil
