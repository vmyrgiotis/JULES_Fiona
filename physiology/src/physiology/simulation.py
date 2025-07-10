import numpy as np
from config import SIM_SETTINGS, PFT_PARAMS
from equations import rate_limiters

def generate_forcings(settings):
    days = settings['days']
    dt = settings['dt_hours']
    t = np.arange(0, days * 24, dt)

    T = settings['T_mean'] + settings['T_amp'] * np.sin(2 * np.pi * t / 24)
    I_par = np.maximum(0, settings['I_max'] * np.sin(np.pi * t / 24))

    P = settings['P']
    ca = settings['ca_ppm'] * 1e-6 * P
    ci = np.full_like(t, settings['chi'] * ca)
    O2 = np.full_like(t, settings['O2_fraction'] * P)
    return t, T, I_par, ci, O2


def run_simulation(pft_key):
    p = PFT_PARAMS[pft_key]
    # Compute Vcmax25 from n0 and n_e, Equation 6
    p['Vcmax25'] = p['n0'] * p['n_e']

    # Generate forcing time series
    t, T, I_par, ci, O2 = generate_forcings(SIM_SETTINGS)

    # Preallocate output
    results = {k: np.zeros_like(t) for k in ['Wc','Wl','We','Rd','Wg','Ap']}

    # Main loop
    for i in range(len(t)):
        out = rate_limiters(T[i], I_par[i], ci[i], O2[i], p)
        for k in results:
            results[k][i] = out[k]

    return t, results
