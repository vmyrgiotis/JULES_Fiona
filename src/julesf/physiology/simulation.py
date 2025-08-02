#simulation.py

import numpy as np

from julesf.physiology.parameters import SIM_SETTINGS, PHOTOSYNTHESIS_PFT_PARAMS, CANOPY_PARAMS, RESPIRATION_PARAMS
from julesf.physiology.equations import rate_limiters, big_leaf_photosynthesis, big_leaf_dark_resp, compute_N_pools, maintenance_resp, growth_resp
    
def generate_forcings(settings):
    days = settings['days']
    dt   = settings['dt_hours']
    t    = np.arange(0, days*24, dt)
    t_day = t % 24

    T     = settings['T_mean'] + settings['T_amp'] * np.sin(2*np.pi*t_day/24)
    I_par = np.maximum(0, settings['I_max'] * np.sin(2*np.pi*t_day/24)) * 1e-6

    ca    = settings['ca_ppm'] * 1e-6 * settings['P']
    ci    = np.full_like(t, settings['chi'] * ca)
    O2    = np.full_like(t, settings['O2_fraction'] * settings['P'])
    return t, T, I_par, ci, O2

def run_photosynthesis(pft_key):

    # Compute Vcmax25 from n0 and n_e, Equation 6
    p = PHOTOSYNTHESIS_PFT_PARAMS[pft_key]
    p['Vcmax25'] = p['n0'] * p['n_e']

    # Generate forcing time series
    t, T, I_par, ci, O2 = generate_forcings(settings=SIM_SETTINGS)

    # Preallocate output
    results = {k: np.zeros_like(t) for k in ['Wc','Wl','We','Rd','Wg','Ap']}

    # Main loop
    for i in range(len(t)):
        out = rate_limiters(T[i], I_par[i], ci[i], O2[i], p)
        for k in results:
            results[k][i] = out[k]

    # big-leaf scaling from leaf to canopy
    k_ext = CANOPY_PARAMS['k_ext']
    LAI   = CANOPY_PARAMS['LAI']
    results['Ac'] = big_leaf_photosynthesis(results['Ap'], k_ext, LAI)

    return t, results


def run_respiration(pft_key):
    """
    EQquation 38, 40-45
    From leaf & canopy photosynthesis + leaf respiration, compute:
    R_dc (canopy dark resp),
    R_pm (maintenance resp),
    R_pg (growth resp),
    R_p  (total plant resp),
    Pi_G (gross primary productivity).
    """
    t, ps = run_photosynthesis(pft_key)

    # canopy dark respiration, Equation 19
    k_ext = CANOPY_PARAMS['k_ext']
    LAI   = CANOPY_PARAMS['LAI']
    R_dc = big_leaf_dark_resp(ps['Rd'], k_ext, LAI)

    # N‐pool ratio, Equation 43-35
    _, _, _, nr_ns_over_nl = compute_N_pools(
        LAI,
        RESPIRATION_PARAMS['n_m'],
        RESPIRATION_PARAMS['sigma1'],
        RESPIRATION_PARAMS['mu_root_leaf_N'],
        RESPIRATION_PARAMS['mu_stem_leaf_N'],
        RESPIRATION_PARAMS['eta_root_C'],
        RESPIRATION_PARAMS['eta_stem_C'],
        CANOPY_PARAMS['h']
    )

    beta = RESPIRATION_PARAMS['beta']
    rg   = RESPIRATION_PARAMS['rg']

    # gross primary productivity, Equation 38
    Pi_G = ps['Ac'] + beta * R_dc

    # maintenance & growth respiration, Equation 40-42
    R_pm = maintenance_resp(R_dc, beta, nr_ns_over_nl)
    R_pg = growth_resp(Pi_G, R_pm, rg)
    R_p  = R_pm + R_pg

    resp = {
        'R_dc':  R_dc,
        'Pi_G':  Pi_G,
        'R_pm':  R_pm,
        'R_pg':  R_pg,
        'R_p':   R_p
    }

    return t, resp

def run_npp(pft_key):
    """
    Equation 39, net primary productivity Pi = Pi_G - R_p.
    """
    t, resp = run_respiration(pft_key)

    Pi_net = resp['Pi_G'] - resp['R_p']
    # bundle everything
    npp = {**resp, 'Pi_net': Pi_net}

    return t, npp

