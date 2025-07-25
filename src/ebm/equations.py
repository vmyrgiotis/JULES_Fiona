import numpy as np

try:
    from parameters import (sigma, rho, c_p, L_v,
                        alpha, epsilon_s,
                        r_a, r_s, r_acap,
                        lambda_soil, C_surf, C_soil1)
except (ModuleNotFoundError):
    from .parameters import (sigma, rho, c_p, L_v,
                        alpha, epsilon_s,
                        r_a, r_s, r_acap,
                        lambda_soil, C_surf, C_soil1)

def H_flux(T_surf, T_air):
    """Sensible heat flux (Eq. 2)."""
    return (rho * c_p / r_a) * (T_surf - T_air)

def E_flux(T_surf, Q1):
    """Latent heat flux (Eq. 3), in energy units W/m2."""
    return (rho / (r_a + r_s)) * (qsat(T_surf) - Q1) * L_v

def qsat(T):
    """Approximate saturation specific humidity (kg/kg) at temp T (K)."""
    return 0.01 * np.exp(0.07 * (T - 273.15))

def G_flux(T_surf, T_soil, nu):
    """Canopy/soil exchange (Eq. 4)."""
    # radiative term
    rad_term  = sigma * epsilon_s * (T_surf**4 - T_soil**4)
    # turbulent term
    turb_term = (rho * c_p / r_acap) * (T_surf - T_soil)
    # soil conduction
    cond_term = lambda_soil * (T_surf - T_soil)
    return nu * (rad_term + turb_term) + (1 - nu) * cond_term

def ebm_rhs(t, T_flat, drivers):
    """
    RHS of surface energy balance, Eq. 1

    T_flat  : [T_surf, ]
    drivers : dict with 'Sdn', 'Ldn', 'Tair', 'Q1', 'nu', 'Tsoil'
    """
    T_surf = T_flat[0]

    Sdn   = drivers["Sdn"](t)
    Ldn   = drivers["Ldn"](t)
    Tair  = drivers["Tair"](t)
    Q1    = drivers["Q1"](t)
    nu    = drivers["nu"](t)
    Tsoil = drivers["Tsoil"](t)

    # net radiation, Eq.1 term 1-3
    Rn = (1 - alpha) * Sdn + epsilon_s * Ldn - sigma * epsilon_s * T_surf**4

    # fluxes, Eq. 1 term 4-6
    H = H_flux(T_surf, Tair)
    E = E_flux(T_surf, Q1)
    G = G_flux(T_surf, Tsoil, nu)

    # surface temperature tendency, Eq. 1
    dTs_dt = (Rn - H - E - G) / C_surf

    return np.array([dTs_dt])
