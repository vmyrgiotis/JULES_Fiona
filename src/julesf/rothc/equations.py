# equations.py

import numpy as np
from src.rothc.parameters import POOLS, kappa, alpha_dr, beta_R, beta_R_default

def F_T_Q10(T_C, Q10=2.0, T_ref=10.0):
    """Temperature modifier 1: Q10 temperature response (Eq.65)."""
    return Q10 ** ((T_C - T_ref) / 10.0)

def F_T_RothC(T_C):
    """Temperature modifier 2: RothC response from Jenkins (1990) (Eq. 66)."""
    return 47.9 / (1.0 + np.exp(106.0 / (T_C + 18.3)))

def F_S(s, s_min=0.1, s_o=0.5):
    """Moisture modifier (Eq. 67)."""
    if s > s_o:
        return 1.0 - 0.8 * (s_o - s) / (1.0 - s_o)
    elif s_min <= s <= s_o:
        return 0.2 + 0.8 * (s - s_min) / (s_o - s_min)
    else:
        return 0.2

def F_V(nu):
    """Vegetation‐cover modifier (Eq. 68)."""
    return 0.6 + 0.4 * (1.0 - nu)

def soil_carbon_rhs(t, C_flat, drivers, texture="loam", temp_fun="Q10", pft="tree"):
    """
    RHS of multi‐pool soil carbon ODEs (Eq.62, 64, 70-73).
    t         : time (days)
    C_flat    : 'DPM', 'RPM', 'BIO', 'HUM'
    drivers   : 'Lambda_c', 'T_soil', 's', 'nu'
    texture   : 'sand', 'loam', 'clay', or 'peat'
    temp_fun  : 'Q10' or 'RothC'
    pft       : 'tree','shrub','grass','crop'

    Returns dC/dt, shape=(4,).
    """
    C = np.asarray(C_flat)
    Λc = drivers["Lambda_c"]
    T  = drivers["T_soil"]
    s  = drivers["s"]
    nu = drivers["nu"]


    # respiration rates, Eq. 64
    FT = F_T_Q10(T) if temp_fun=="Q10" else F_T_RothC(T)
    FS = F_S(s)
    FV = F_V(nu)
    R_i = np.array([kappa[p] * C[i] * FT * FS * FV
                    for i, p in enumerate(POOLS)])
    R_s = R_i.sum()

    # litter partitioning, Eq. 69
    α     = alpha_dr[pft]
    f_dpm = α / (1 + α)
    f_rpm = 1 - f_dpm

    # build dC/dt according to Eq.70-733
    β = beta_R.get(texture, beta_R_default)
    dC = np.zeros_like(C)
    dC[0] = f_dpm * Λc       - R_i[0]       # DPM
    dC[1] = f_rpm * Λc       - R_i[1]       # RPM
    dC[2] = 0.46 * β * R_s   - R_i[2]       # BIO
    dC[3] = 0.54 * β * R_s   - R_i[3]       # HUM

    return dC
