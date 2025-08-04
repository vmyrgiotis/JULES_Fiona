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
    if soil_initial and 'T_soil' in soil_initial:
        mean_tsoil = np.mean(soil_initial['T_soil'])
        ebm_drivers_complete['Tsoil'] = lambda t: mean_tsoil
    else:
        ebm_drivers_complete['Tsoil'] = lambda t: ebm_drivers_complete['Tair'](t) - 2.0
    
    print(f"  ✓ EBM coupling: nu={nu_cover:.3f}")
    
    # ========== 2. RUN EBM ==========
    T0 = soil_initial['T_soil'][0] if (soil_initial and 'T_soil' in soil_initial) else ebm_drivers_complete['Tair'](0)
    print(f"  ✓ EBM initial temperature: {T0:.1f}K")
    
    t_ebm, Ts_ebm = run_ebm(t_span, T0, ebm_drivers_complete, dt_out=dt_hours/24)
    print(f"  ✓ EBM complete: T_surface range {min(Ts_ebm):.1f}-{max(Ts_ebm):.1f}K")
    
    # ========== 3. SETUP PHYSIOLOGY DRIVERS ==========
    physiology_forcing = {}
    
    if external_drivers and 'era5_forcing' in external_drivers:
        era5 = external_drivers['era5_forcing']
        physiology_forcing['pressure'] = lambda t: era5['pressure'](t)  # Pa
        physiology_forcing['co2'] = lambda t: era5['co2'](t)  # ppm
        print(f"  ✓ Physiology: Using ERA5 pressure and CO2")
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
        soil_drivers_complete['precip'] = lambda t: max(0, era5['precip'](t))  # kg/m²/s
        
        # Test precipitation
        try:
            test_precip = soil_drivers_complete['precip'](24.0)  # 1 day in
            print(f"  ✓ Soil: Using ERA5 precipitation (test: {test_precip:.6f} kg/m²/s)")
        except Exception as e:
            print(f"  ❌ ERA5 precipitation test failed: {e}")
    else:
        soil_drivers_complete['precip'] = lambda t: 1e-5  # Default precipitation
        print(f"  ✓ Soil: Using default precipitation")
    
    # Add EBM surface temperature coupling
    soil_drivers_complete['T_surface'] = lambda t: np.interp(t, t_ebm, Ts_ebm)
    
    # Add vegetation coupling
    soil_drivers_complete['nu'] = lambda t: nu_cover  # For evapotranspiration
    
    # ========== 6. RUN SOIL MODEL ==========
    t_soil, theta_ts, T_ts = solve_soil_rhs(
        t_span=(0, days*24),
        theta_init=soil_initial['theta'] if soil_initial else None,
        T_init=soil_initial['T_soil'] if soil_initial else None,
        drivers=soil_drivers_complete,
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
    """
    Run physiology with coupled inputs from EBM, TRIFFID, and Soil
    """
    from julesf.physiology.parameters import SIM_SETTINGS
    from julesf.physiology.simulation import generate_forcings
    
    def gen_forcings_coupled(settings):
        t = np.arange(0, days * 24, dt_hours)
        
        # Temperature from EBM (coupled)
        T_leaf_C = np.array([physiology_forcing['T_leaf'](ti) for ti in t])
        
        # I_par: SMART BACKUP PLAN - simple sinusoidal pattern
        # TODO: Implement JULES original method
        hour_of_day = t % 24
        daylight_hours = (hour_of_day >= 6) & (hour_of_day <= 18)
        I_par = np.where(daylight_hours, 
                        800 * np.sin(np.pi * (hour_of_day - 6) / 12),  # Peak at noon
                        0)
        I_par = np.maximum(I_par, 0)  # No negative light
        
        # CO2 from ERA5 (external forcing)
        ci = np.array([physiology_forcing['co2'](ti) for ti in t])
        
        # Oxygen (constant)
        O2 = np.full_like(t, 21000)  # ppmv
        
        return t, T_leaf_C, I_par, ci, O2
    
    # Monkey patch and run
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