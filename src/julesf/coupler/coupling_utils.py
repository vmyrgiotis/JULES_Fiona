"""
Utilities for coupling JULES component models
"""
import numpy as np

def convert_npp_to_carbon(npp_umol_m2_s, integration_hours=168):
    """Convert μmol CO₂ m⁻² s⁻¹ (from physiology model) to kg C m⁻² week⁻¹ (for TRIFFID)"""
    mol_co2_per_s = npp_umol_m2_s / 1e6
    kg_c_per_s = mol_co2_per_s * 0.012
    kg_c_per_period = kg_c_per_s * 3600 * integration_hours
    return kg_c_per_period

def extract_vegetation_vars(triffid_state):
    """Extract coupling variables from TRIFFID state vector"""
    Lb_tree, Lb_grass, nu_tree, nu_grass = triffid_state
    
    return {
        'LAI_total': Lb_tree + Lb_grass,
        'nu_total': nu_tree + nu_grass,
        'Lb_tree': Lb_tree,
        'Lb_grass': Lb_grass,
        'nu_tree': nu_tree,
        'nu_grass': nu_grass
    }

def weekly_mean_soil_moisture(soil_theta_timeseries):
    """Compute weekly mean soil moisture for RothC"""
    return np.mean(soil_theta_timeseries[0, :])  # Surface layer

def calculate_litterfall(triffid_state, pft_params):
    """Calculate litterfall from TRIFFID state (kg C m⁻² day⁻¹) (for coupling to RothC)"""
    Lb_tree, Lb_grass = triffid_state[0], triffid_state[1]
    
    sigma1_tree = pft_params['broadleaf_tree']['sigma1']
    sigma1_grass = pft_params['C3_grass']['sigma1']
    gamma_l_tree = pft_params['broadleaf_tree']['gamma_l']
    gamma_l_grass = pft_params['C3_grass']['gamma_l']
    
    L_tree = sigma1_tree * Lb_tree
    L_grass = sigma1_grass * Lb_grass
    
    # Annual rate converted to daily
    litterfall = (gamma_l_tree * L_tree + gamma_l_grass * L_grass) / 365.0
    
    return litterfall