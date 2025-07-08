import numpy as np
from equation import (
    compute_transmissivity, compute_baseflow,
    compute_lambda_c, compute_f_sat, update_water_table
)

def simulate_topmodel(
    nx=5, ny=5, days=365, dt=24*3600,
    T0=1e-5, f=0.2, a_s=0.3, c_s=1.0, theta_eff=0.2
):
    # set up grid & parameters
    np.random.seed(0)
    lambda_bar = np.random.uniform(3.0, 6.0, (nx, ny))
    Rb_max     = T0 * np.exp(-lambda_bar)

    # forcing - a sine wave with peak during mid-year
    t   = np.arange(days)
    W0m = 1.0 / (365 * 86400)
    W0  = W0m * (1 + 0.5 * np.sin(2*np.pi*t/days))

    # initial state & diagnostics
    z_bar    = np.full((days, nx, ny), 2.0)
    Rb       = np.zeros_like(z_bar)
    lambda_c = np.zeros_like(z_bar)
    f_sat    = np.zeros_like(z_bar)
    Rse      = np.zeros_like(z_bar)

    # main time loop
    for i in range(days-1):
        Tz        = compute_transmissivity(z_bar[i], T0, f)
        Rb[i]     = compute_baseflow(Tz, lambda_bar)
        lambda_c[i] = compute_lambda_c(Rb[i], Rb_max)
        f_sat[i]  = compute_f_sat(lambda_c[i], a_s, c_s)
        Rse[i]    = f_sat[i] * W0[i]
        z_bar[i+1]= update_water_table(z_bar[i], W0[i], Rb[i], dt, theta_eff)

    # final diagnostics
    Tz        = compute_transmissivity(z_bar[-1], T0, f)
    Rb[-1]    = compute_baseflow(Tz, lambda_bar)
    lambda_c[-1] = compute_lambda_c(Rb[-1], Rb_max)
    f_sat[-1] = compute_f_sat(lambda_c[-1], a_s, c_s)
    Rse[-1]   = f_sat[-1] * W0[-1]

    return t, lambda_bar, W0, z_bar, Rb, lambda_c, f_sat, Rse
