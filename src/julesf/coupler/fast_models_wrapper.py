"""
Wrapper for fast models (EBM + Soil + Physiology) at 0.5h timestep
"""
import numpy as np
from scipy.integrate import solve_ivp

from julesf.ebm.run_ebm import run_ebm, drivers as ebm_drivers
from julesf.physiology.simulation import run_npp
from julesf.soil.simulation import solve_soil_rhs
from julesf.coupler.coupling_utils import convert_npp_to_carbon

def run_fast_models_week(nu_cover, LAI_total, soil_initial, week_num=0):
    """Run EBM + Soil + Physiology for 1 week at 0.5h timestep"""
    days = 7
    dt_hours = 0.5
    t_span = (0, days)
    
    print(f"Week {week_num+1}: Running fast models (EBM+Soil+Physiology)...")
    
    # 1. Setup EBM with soil feedbacks
    ebm_drivers_mod = ebm_drivers.copy()
    ebm_drivers_mod['nu'] = lambda t: nu_cover
    
    # 2. Run EBM 
    T0 = soil_initial['T_soil'][0] if soil_initial else ebm_drivers['Tair'](0)
    t_ebm, Ts_ebm = run_ebm(t_span, T0, ebm_drivers_mod, dt_out=dt_hours/24)
    
    # 3. Run Physiology with EBM temperature
    t_phys, npp = run_physiology_with_ebm(Ts_ebm, LAI_total, days, dt_hours)
    
    # 4. Run Soil model with EBM surface temperature
    soil_drivers = {
        'T_surface': lambda t: np.interp(t, t_ebm, Ts_ebm),
        'precip': lambda t: 1e-5,  # kg/m²/s, placeholder
        'ET': lambda t: 5e-6,      # kg/m²/s, placeholder
    }
    
    t_soil, theta_ts, T_ts = solve_soil_rhs(
        t_span=(0, days*24),
        theta_init=soil_initial['theta'] if soil_initial else None,
        T_init=soil_initial['T_soil'] if soil_initial else None,
        drivers=soil_drivers,
        dt_out=dt_hours
    )
    print(f"Final moisture range: {theta_ts.min():.3f} - {theta_ts.max():.3f} m³/m³")
    print(f"Final temperature range: {T_ts.min():.1f} - {T_ts.max():.1f} K")
    
    soil_results = {
        'theta': theta_ts,                # n_layers × n_times
        'T_soil': T_ts                    # n_layers × n_times
    }
    
    # 5. Compute integrated NPP
    npp_mean = np.mean(npp['Pi_net'])  # μmol CO₂ m⁻² s⁻¹
    npp_for_triffid = convert_npp_to_carbon(npp_mean, days * 24)

    Pi_G_mean = np.mean(npp.get('Pi_G', npp['Pi_net']))  # Gross photosynthesis
    R_p_mean = np.mean(npp.get('R_p', 0))                # Plant respiration
    
    
    # 6. Return results
    return {
        'T_surface': Ts_ebm,                  # K, time series
        'soil_theta': soil_results['theta'],  # m³/m³, time series  
        'soil_T': soil_results['T_soil'],     # K, time series
        'NPP_total': npp_for_triffid,         # kg C m⁻² week⁻¹
        'Pi_G': Pi_G_mean,                    # μmol CO₂ m⁻² s⁻¹
        'R_p': R_p_mean,                      # μmol CO₂ m⁻² s⁻¹
        'soil_final': {                       # Final states for next week
            'theta': soil_results['theta'][:, -1],
            'T_soil': soil_results['T_soil'][:, -1]
        }
    }

def run_physiology_with_ebm(Ts_ebm, LAI, days=7, dt_hours=0.5):
    """Run physiology with EBM surface temperature"""
    from julesf.physiology.parameters import SIM_SETTINGS
    from julesf.physiology.simulation import generate_forcings
    
    # Override original forcing generator to use EBM temperature
    def gen_forcings_override(settings):
        t = np.arange(0, days * 24, dt_hours)
        T_leaf_C = Ts_ebm - 273.15  # K to °C
        _, _, I_par, ci, O2 = generate_forcings(settings)
        return t, T_leaf_C, I_par, ci, O2
    
    # Monkey patch the forcing generator
    import julesf.physiology.simulation as sim
    original_gen_forcings = sim.generate_forcings
    sim.generate_forcings = gen_forcings_override
    
    # Run NPP calculation with modified settings
    settings_copy = SIM_SETTINGS.copy()
    settings_copy['days'] = days
    settings_copy['dt_hours'] = dt_hours
    
    # Use needleleaf_tree as default PFT (should be configurable)
    t_npp, npp = run_npp('needleleaf_tree')
    
    # Restore original function
    sim.generate_forcings = original_gen_forcings
    
    return t_npp, npp