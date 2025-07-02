import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# parameters for two PFTs: broadleaf tree and C3 grass, from Clark et al (2011) Table 6 and 7
params = {
    'sigma1': np.array([0.0375, 0.0250]),  # leaf C per LAI
    'awl': np.array([0.65, 0.005]),        # woody allometry coeff
    'bwl': np.array([1.667, 1.667]),       # woody allometry exponent
    'gamma_l': np.array([0.25, 0.25]),     # leaf turnover (yr^-1)
    'gamma_r': np.array([0.25, 0.25]),     # root turnover
    'gamma_w': np.array([0.005, 0.20]),    # wood turnover
    'Lmin': np.array([1.0, 1.0]),          # min LAI
    'Lmax': np.array([9.0, 4.0]),          # max LAI
    'c': np.array([[1.0, 1.0], [0.0, 1.0]]),
    # competition coefficient - the first matrix for intra-species competition, the second for inter-species competition
}


# partitioning function, according to Eq. 54
def lambda_fun(Lb, Lmin, Lmax):
    lam = np.clip((Lb - Lmin) / (Lmax - Lmin), 0, 1)
    return lam


# generating a sine wave as my artificial NPP forcing
def Pi_t(t, A=0.01, B=0.01, period=10.0):
    return A * np.sin(2 * np.pi * t / period) + B


# right-hand-side of ODEs
def triffid_rhs(t, y):
    Lb = y[:2]  # first variable we are tracking, balanced LAI (of both PFT)
    nu = y[2:]  # second variable, functional cover of both PFT

    # leaf, root, and wood biomass pools, according to Eq. 56-58
    L = params['sigma1'] * Lb
    R = L
    W = params['awl'] * Lb ** params['bwl']
    Cv = L + R + W

    # lambda_i, as defined above
    lam = lambda_fun(Lb, params['Lmin'], params['Lmax'])

    # NPP for each PFT, using the pre-defined sine wave
    Pi = np.array([Pi_t(t), Pi_t(t)])

    # dLb/dt, derived from Eq. 51 (This one is a bit tricky - I am still working on the maths!)
    dLb_dt = (1 - lam) * Pi / params['sigma1'] \
             - (params['gamma_l'] + params['gamma_r']) * Lb \
             - (params['gamma_w'] * W) / params['sigma1']

    # dnu/dt, derived from Eq. 52
    competition = 1 - (params['c'] @ nu)
    dnu_dt = lam * Pi * nu * competition / Cv

    return np.concatenate([dLb_dt, dnu_dt])

# initial condition
y0 = [9.0, 3.0, 0.5, 0.5]
t_span = (0, 50)
t_eval = np.linspace(0, 50, 501)

# integration
sol = solve_ivp(triffid_rhs, t_span, y0, t_eval=t_eval, method='RK45')

# visualization
fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
axes[0].plot(sol.t, sol.y[2], label='Tree cover')
axes[0].plot(sol.t, sol.y[3], label='Grass cover')
axes[0].set_ylabel('Fractional cover')
axes[0].legend()

axes[1].plot(sol.t, sol.y[0], label='Tree Lb')
axes[1].plot(sol.t, sol.y[1], label='Grass Lb')
axes[1].set_ylabel('Balanced LAI')
axes[1].set_xlabel('Time (yr)')
axes[1].legend()

plt.tight_layout()
plt.show()
