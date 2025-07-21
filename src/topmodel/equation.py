import numpy as np

def compute_transmissivity(z, T0, f):
    return T0 * np.exp(-f * z)

def compute_baseflow(Tz, lambda_bar):
    return Tz * np.exp(-lambda_bar)

def compute_lambda_c(Rb, Rb_max):
    return np.log(Rb_max / Rb)

def compute_f_sat(lambda_c, a_s, c_s):
    return a_s * np.exp(-c_s * lambda_c)

def update_water_table(z, W0, Rb, dt, theta_eff):
    dz = (W0 - Rb) * dt / theta_eff
    return np.maximum(z + dz, 0.0)
