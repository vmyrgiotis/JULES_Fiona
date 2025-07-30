# equations.py (Equations 51-59 from Clark et al 2011)
import numpy as np

def biomass_pools(Lb, params):
    """Calculate leaf, root, wood biomass from balanced LAI (Eq. 56-58)"""
    L = params['sigma1'] * Lb  # kg C m⁻²
    R = L                      # kg C m⁻² (root = leaf)
    W = params['awl'] * Lb ** params['bwl']  # kg C m⁻²
    Cv = L + R + W            # kg C m⁻², total vegetation carbon
    return L, R, W, Cv

def litterfall_rate(L, R, W, params):
    """Local litterfall rate (Eq. 59)"""
    Lambda1 = (params['gamma_l'] * L + 
               params['gamma_r'] * R + 
               params['gamma_w'] * W)  # kg C m⁻² yr⁻¹
    return Lambda1

def partitioning_function(Lb, Lmin, Lmax):
    """Carbon partitioning function (Eq. 54)"""
    lam = np.clip((Lb - Lmin) / (Lmax - Lmin), 0, 1)
    return lam

def triffid_rhs(t, y, params_combined, npp_external=None):
    """
    TRIFFID right-hand side (Equations 51-52)
    
    Args:
        t: time (years)
        y: state vector [Lb_tree, Lb_grass, nu_tree, nu_grass]
        params_combined: Combined parameters dict
        npp_external: External NPP [kg C m⁻² yr⁻¹] or None for artificial
    """
    Lb = y[:2]  # Balanced LAI [m² m⁻²]
    nu = y[2:]  # Fractional cover [dimensionless]
    
    # Extract PFT parameters
    sigma1 = np.array([params_combined['broadleaf_tree']['sigma1'],
                       params_combined['C3_grass']['sigma1']])
    awl = np.array([params_combined['broadleaf_tree']['awl'],
                    params_combined['C3_grass']['awl']])
    bwl = np.array([params_combined['broadleaf_tree']['bwl'],
                    params_combined['C3_grass']['bwl']])
    Lmin = np.array([params_combined['broadleaf_tree']['Lmin'],
                     params_combined['C3_grass']['Lmin']])
    Lmax = np.array([params_combined['broadleaf_tree']['Lmax'],
                     params_combined['C3_grass']['Lmax']])
    
    # Biomass pools
    L, R, W, Cv = biomass_pools(Lb, {'sigma1': sigma1, 'awl': awl, 'bwl': bwl})
    
    # Litterfall
    Lambda1 = litterfall_rate(L, R, W, {
        'gamma_l': np.array([params_combined['broadleaf_tree']['gamma_l'],
                             params_combined['C3_grass']['gamma_l']]),
        'gamma_r': np.array([params_combined['broadleaf_tree']['gamma_r'],
                             params_combined['C3_grass']['gamma_r']]),
        'gamma_w': np.array([params_combined['broadleaf_tree']['gamma_w'],
                             params_combined['C3_grass']['gamma_w']])
    })
    
    # Partitioning function
    lam = partitioning_function(Lb, Lmin, Lmax)
    
    # NPP input
    if npp_external is not None:
        Pi = npp_external  # kg C m⁻² yr⁻¹
    else:
        # Artificial NPP for standalone testing
        Pi = np.array([artificial_npp(t), artificial_npp(t)])
    
    # dLb/dt (Equation 51)
    numerator = (1 - lam) * Pi - Lambda1
    denominator = 2 * sigma1 + awl * bwl * Lb**(bwl - 1)
    dLb_dt = numerator / denominator
    
    # dnu/dt (Equation 52)
    c_matrix = params_combined['competition']['c']
    competition = 1 - (c_matrix @ nu)
    dnu_dt = lam * Pi * nu * competition / Cv
    
    return np.concatenate([dLb_dt, dnu_dt])

def artificial_npp(t, A=0.1, B=0.1, period=5.0):
    """Artificial NPP for standalone testing"""
    return A * np.sin(2 * np.pi * t / period) + B  # kg C m⁻² yr⁻¹