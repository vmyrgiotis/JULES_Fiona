#equations.py
import numpy as np

def temperature_response(T, Vcmax25, Q10_Kc, T_low, T_upp):
    """
    Equation 4, 5
    """
    f_T = Q10_Kc ** ((T - 25.0) / 10.0)
    denom = (1.0 + np.exp(0.3 * (T - T_upp))) * (1.0 + np.exp(0.3 * (T_low - T)))
    return Vcmax25 * f_T / denom


def michaelis_constants(T, Q10_Kc, Q10_Ko):
    """
    Equation 9, 10
    """
    Kc = 30.0 * Q10_Kc ** ((T - 25.0) / 10.0)
    Ko = 3e4  * Q10_Ko ** ((T - 25.0) / 10.0)
    return Kc, Ko


def compensation_point(T, O2, Q10_rs):
    """
    Γ, the CO2 compensation point, Equation 7, 8
    """
    tau = 2600.0 * Q10_rs ** ((T - 25.0) / 10.0)
    return O2 / (2.0 * tau)


def rate_limiters(T, I_par, ci, O2, params):
    """
    Compute Wc, Wl, We, according to Equation 1-3, and 13-14
    alpha, omega, f_dr, Vcmax25, Q10_Kc, Q10_Ko, Q10_rs, T_low, T_upp
    """

    alpha = params['alpha']
    omega = params['omega']
    f_dr  = params['f_dr']
    Vcmax25 = params['Vcmax25']
    Q10_Kc, Q10_Ko, Q10_rs = params['Q10_Kc'], params['Q10_Ko'], params['Q10_rs']
    T_low, T_upp = params['T_low'], params['T_upp']

    # Compute Vcmax
    Vcmax = temperature_response(T, Vcmax25, Q10_Kc, T_low, T_upp)
    # Kc, Ko, Gamma
    Kc, Ko = michaelis_constants(T, Q10_Kc, Q10_Ko)
    Gamma = compensation_point(T, O2, Q10_rs)

    # Limiters
    Wc = Vcmax * (ci - Gamma) / (ci + Kc * (1.0 + O2 / Ko))
    Wl = alpha * (1.0 - omega) * I_par * (ci - Gamma) / (ci + 2.0 * Gamma)
    We = 0.5 * Vcmax

    # Respiration and net potential photosynthesis
    Rd = f_dr * Vcmax                     # Equation 13
    Wg = np.minimum.reduce([Wc, Wl, We])  # simplified from equation 11 and 12
    Ap = Wg - Rd                          # Equation 14


    return dict(Wc=Wc, Wl=Wl, We=We, Rd=Rd, Wg=Wg, Ap=Ap)

def big_leaf_photosynthesis(A_leaf, k_ext, LAI):
    """
    Equation 19. The second approach for leaf-to-canopy scaling will be updated later.
    """
    if k_ext <= 0:
        return A_leaf * LAI
    return A_leaf * (1.0 - np.exp(-k_ext * LAI)) / k_ext


def compute_N_pools(LAI, n_m, sigma1,
                    mu_rl, mu_sl,
                    eta_root_C, eta_stem_C, h):
    """
    Equation 43-45
    N_l = n_m * sigma1 * LAI
    N_r = mu_rl * n_m * (eta_root_C * LAI)
    N_s = mu_sl * n_m * (eta_stem_C * h * LAI)
    """
    N_l = n_m * sigma1 * LAI
    N_r = mu_rl * n_m * (eta_root_C * LAI)
    N_s = mu_sl * n_m * (eta_stem_C * h * LAI)
    ratio = (N_r + N_s) / N_l if N_l > 0 else 0.0
    return N_l, N_r, N_s, ratio


def big_leaf_dark_resp(Rd_leaf, k_exT, LAI):
    """Scale leaf Rd to canopy dark respiration R_dc."""
    return big_leaf_photosynthesis(Rd_leaf, k_exT, LAI)


def maintenance_resp(R_dc, beta, nr_ns_over_nl):
    """Equation 42, R_pm = 0.012 * R_dc * (beta + (N_r+N_s)/N_l)"""
    return 0.012 * R_dc * (beta + nr_ns_over_nl)


def growth_resp(Pi_G, R_pm, rg):
    """Equation 41, R_pg = rg * (Pi_G - R_pm)"""
    return rg * (Pi_G - R_pm)
