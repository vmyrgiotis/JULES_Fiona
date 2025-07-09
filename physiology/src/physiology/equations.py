import numpy as np

def compute_temperature_response(T, Vcmax25, Q10_Kc, T_low, T_upp):
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
    Vcmax = compute_temperature_response(T, Vcmax25, Q10_Kc, T_low, T_upp)
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