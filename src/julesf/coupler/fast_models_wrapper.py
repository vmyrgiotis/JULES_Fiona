"""
Wrapper for fast models (EBM + Soil + Physiology) at 0.5h timestep
"""
import numpy as np
from scipy.integrate import solve_ivp

from julesf.ebm.run_ebm import run_ebm
from julesf.physiology.simulation import run_npp
from julesf.soil.simulation import solve_soil_rhs
from julesf.coupler.coupling_utils import convert_npp_to_carbon

def run_fast_models_week(nu_cover, LAI_total, soil_initial, week_num=0, external_drivers=None):
    """
    Run EBM + Soil + Physiology for 1 week at 0.5h timestep
    
    CLEAR INPUT STRUCTURE:
    - external_drivers: Optional external forcing {'era5_forcing': {...}}
    - nu_cover, LAI_total: from TRIFFID coupling
    - soil_initial: from previous week's soil state
    
    EXTERNAL vs INTERNAL INPUTS:
    - EBM: ERA5(Tair,Sdn,Ldn,Q1) + COUPLING(nu,Tsoil) 
    - Physiology: ERA5(pressure,co2) + EBM(T_leaf) + TRIFFID(LAI,nu) + SOIL(beta)
    - Soil: ERA5(precip) + EBM(T_surface) + TRIFFID(nu)
    """
    days = 7
    dt_hours = 0.5
    t_span = (0, days)
    
    print(f"Week {week_num+1}: Running fast models (EBM+Soil+Physiology)...")
    
    # ========== DEBUG: CHECK INPUTS ==========
    if external_drivers:
        print(f"  External drivers: {list(external_drivers.keys())}")
        if 'era5_forcing' in external_drivers:
            era5_vars = list(external_drivers['era5_forcing'].keys())
            print(f"  ERA5 variables: {era5_vars}")
        else:
            print(f"  ⚠️  ERROR: 'era5_forcing' not found in external_drivers!")
    else:
        print(f"  Using internal synthetic forcing")
    
    # ========== 1. SETUP EBM DRIVERS ==========
    ebm_drivers_complete = {}
    
    if external_drivers and 'era5_forcing' in external_drivers:
        # ERA5 meteorological forcing
        era5 = external_drivers['era5_forcing']
        
        # Required ERA5 variables for EBM
        required_vars = ['Tair', 'Sdn', 'Ldn', 'Q1']
        missing_vars = [var for var in required_vars if var not in era5]
        
        if missing_vars:
            print(f"  ⚠️  WARNING: Missing ERA5 variables for EBM: {missing_vars}")
        
        ebm_drivers_complete['Tair'] = lambda t: era5['Tair'](t * 24)  # Convert days to hours
        ebm_drivers_complete['Sdn'] = lambda t: max(0, era5['Sdn'](t * 24))
        ebm_drivers_complete['Ldn'] = lambda t: era5['Ldn'](t * 24)
        ebm_drivers_complete['Q1'] = lambda t: era5['Q1'](t * 24)  # Estimated humidity
        
        # Test ERA5 drivers
        try:
            test_t = 0.1  # 0.1 days = 2.4 hours
            test_temp = ebm_drivers_complete['Tair'](test_t)
            test_sdn = ebm_drivers_complete['Sdn'](test_t)
            print(f"  ✓ ERA5 test: Tair={test_temp:.1f}K, Sdn={test_sdn:.1f}W/m²")
        except Exception as e:
            print(f"  ❌ ERA5 driver test failed: {e}")
            
        print(f"  ✓ EBM: Using ERA5 meteorological forcing")
        
    else:
        # Use default synthetic drivers for ERA5 variables
        from julesf.ebm.run_ebm import drivers as default_ebm
        ebm_drivers_complete['Tair'] = default_ebm['Tair']
        ebm_drivers_complete['Sdn'] = default_ebm['Sdn']
        ebm_drivers_complete['Ldn'] = default_ebm['Ldn']
        ebm_drivers_complete['Q1'] = default_ebm['Q1']
        
        # Test synthetic drivers
        try:
            test_temp = ebm_drivers_complete['Tair'](0.1)
            test_sdn = ebm_drivers_complete['Sdn'](0.1)
            print(f"  ✓ Synthetic test: Tair={test_temp:.1f}K, Sdn={test_sdn:.1f}W/m²")
        except Exception as e:
            print(f"  ❌ Synthetic driver test failed: {e}")
            
        print(f"  ✓ EBM: Using synthetic meteorological forcing")
    
    # Add coupling variables (internal)
    ebm_drivers_complete['nu'] = lambda t: nu_cover  # From TRIFFID
    
    # Initial Tsoil estimate (will be updated by soil coupling in subsequent weeks)
    if soil_initial is not None and isinstance(soil_initial, dict) and 'T_soil' in soil_initial:
        mean_tsoil = np.mean(soil_initial['T_soil'])
        ebm_drivers_complete['Tsoil'] = lambda t: mean_tsoil
    else:
        ebm_drivers_complete['Tsoil'] = lambda t: ebm_drivers_complete['Tair'](t) - 2.0

    print(f"  ✓ EBM coupling: nu={nu_cover:.3f}")

    # Initial temperature for EBM
    if soil_initial is not None and isinstance(soil_initial, dict) and 'T_soil' in soil_initial:
        T0 = soil_initial['T_soil'][0]
    else:
        T0 = ebm_drivers_complete['Tair'](0)
    
    print(f"  ✓ EBM initial temperature: {T0:.1f}K")
    
    t_ebm, Ts_ebm = run_ebm(t_span, T0, ebm_drivers_complete, dt_out=dt_hours/24)
    print(f"  ✓ EBM complete: T_surface range {min(Ts_ebm):.1f}-{max(Ts_ebm):.1f}K")
    
    # ========== 3. SETUP PHYSIOLOGY DRIVERS ==========
    physiology_forcing = {}
    
    if external_drivers and 'era5_forcing' in external_drivers:
        era5 = external_drivers['era5_forcing']
        physiology_forcing['pressure'] = lambda t: era5['pressure'](t)  # Pa
        physiology_forcing['co2'] = lambda t: era5['co2'](t)  # ppm
        physiology_forcing['PAR'] = lambda t: era5['PAR'](t)  # W/m²
        physiology_forcing['T_leaf'] = lambda t: np.interp(t, t_ebm, Ts_ebm) - 273.15  # K to °C
        print(f"  ✓ Physiology: Using ERA5 pressure, CO2, and PAR")
    else:
        physiology_forcing['pressure'] = lambda t: 101325.0  # Standard pressure
        physiology_forcing['co2'] = lambda t: 400.0  # Default CO2
        print(f"  ✓ Physiology: Using default pressure and CO2")
    
    # Add EBM surface temperature coupling
    physiology_forcing['T_leaf'] = lambda t: np.interp(t, t_ebm, Ts_ebm) - 273.15  # K to °C
    
    # ========== 4. RUN PHYSIOLOGY ==========
    t_phys, npp = run_physiology_coupled(
        physiology_forcing, LAI_total, nu_cover, soil_initial, days, dt_hours
    )
    print(f"  ✓ Physiology complete: NPP={np.mean(npp['Pi_net']):.6f} μmol CO2/m²/s")
    
    # ========== 5. SETUP SOIL DRIVERS ==========
    soil_drivers_complete = {}
    
    if external_drivers and 'era5_forcing' in external_drivers:
        era5 = external_drivers['era5_forcing']
        soil_drivers_complete['precipitation']    = lambda t: max(0, era5['precip'](t))
        soil_drivers_complete['evapotranspiration'] = lambda t: era5['evapotranspiration'](t)
        print("  ✓ Soil: using ERA5 ET forcing")
                
    else:
        soil_drivers_complete['precipitation'] = lambda t: 1e-5
        soil_drivers_complete['evapotranspiration'] = lambda t: 0.0
        print(f"  ✓ Soil: no ET input, set to zero")

    # Add EBM surface temperature coupling
    soil_drivers_complete['T_surface'] = lambda t: np.interp(t, t_ebm, Ts_ebm)
    
    # Add vegetation coupling
    soil_drivers_complete['nu'] = lambda t: nu_cover  # For evapotranspiration
    
    
    # ========== 6. RUN SOIL MODEL ==========
    theta_arg = None
    T_arg = None
    
    # If soil_initial is provided, use it to initialize soil state
    if soil_initial is not None:
        if isinstance(soil_initial, dict):
            # Dictionary format: {'theta': array, 'T_soil': array}
            theta_arg = soil_initial.get('theta')
            T_arg = soil_initial.get('T_soil')
        else:
            # Assume it's just a theta array
            theta_arg = soil_initial
            T_arg = None  # Will use default in solve_soil_rhs
    
    t_soil, theta_ts, T_ts = solve_soil_rhs(
        t_span=(0, days*24),
        theta_init=theta_arg,
        T_init=T_arg,
        drivers=soil_drivers_complete,
        method="BDF",
        dt_out=dt_hours
    )
    print(f"  ✓ Soil complete: θ={np.mean(theta_ts[0,:]):.3f} m³/m³, T={np.mean(T_ts[0,:])-273.15:.1f}°C")
    
    # ========== 7. PROCESS RESULTS ==========
    npp_mean = np.mean(npp['Pi_net'])
    npp_for_triffid = convert_npp_to_carbon(npp_mean, days * 24)
    
    return {
        'T_surface': Ts_ebm,
        'soil_theta': theta_ts,
        'soil_T': T_ts,
        'NPP_total': npp_for_triffid,
        'Pi_G': np.mean(npp.get('Pi_G', npp['Pi_net'])),
        'R_p': np.mean(npp.get('R_p', 0)),
        'soil_final': {
            'theta': theta_ts[:, -1],
            'T_soil': T_ts[:, -1]
        }
    }

def run_physiology_coupled(physiology_forcing, LAI, nu_cover, soil_state, days=7, dt_hours=0.5):
    from julesf.physiology.parameters import SIM_SETTINGS
    from julesf.physiology.simulation import generate_forcings
    
    def gen_forcings_coupled(settings):
        t = np.arange(0, days * 24, dt_hours)
        T_leaf_C  = np.array([physiology_forcing['T_leaf'](ti) for ti in t])        
        I_par     = np.array([physiology_forcing['PAR'](ti)    for ti in t])
        ci        = np.array([physiology_forcing['co2'](ti)    for ti in t])
        O2        = np.full_like(t, 21000)  # ppmv
        
        return t, T_leaf_C, I_par, ci, O2
    
    import julesf.physiology.simulation as sim
    original_gen_forcings = sim.generate_forcings
    sim.generate_forcings = gen_forcings_coupled
    
    settings_copy = SIM_SETTINGS.copy()
    settings_copy['days'] = days
    settings_copy['dt_hours'] = dt_hours
    
    t_npp, npp = run_npp('needleleaf_tree')
    
    # Restore original function
    sim.generate_forcings = original_gen_forcings
    
    return t_npp, npp

# Previous method for calculating ET (not using ERA5 data), uncomment if needed for future use
# def calculate_blaney_criddle_et(temperature_K, latitude_deg=52.0):
#     """
#     Calculate potential evapotranspiration using Blaney-Criddle method
#     TODO: Replace with JULES original method/external forcing data if available
    
#     PET₀ = p(0.457 T_mean + 8.128) [mm/day]
    
#     Parameters:
#     - temperature_K: air temperature in Kelvin
#     - latitude_deg: latitude in degrees (default 52°N for UK)
    
#     Returns:
#     - PET in kg/m²/s
#     """
#     # Convert K to °C
#     T_celsius = temperature_K - 273.15
    
#     # Simplified p factor for mid-latitudes (varies seasonally, but use average)
#     p_factor = 0.488
    
#     # Blaney-Criddle formula (mm/day)
#     PET_mm_day = p_factor * (0.457 * T_celsius + 8.128)
    
#     # Ensure non-negative
#     PET_mm_day = max(0, PET_mm_day)
    
#     # Convert mm/day to kg/m²/s
#     # 1 mm/day = 1 kg/m²/day = 1/(24×3600) kg/m²/s
#     PET_kg_m2_s = PET_mm_day / (24 * 3600)
    
#     return PET_kg_m2_s
