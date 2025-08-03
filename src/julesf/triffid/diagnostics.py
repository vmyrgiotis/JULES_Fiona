import numpy as np
import matplotlib.pyplot as plt

from julesf.triffid.simulation import solve_triffid
from julesf.triffid.parameters import PFT_PARAMS, COMPETITION_MATRIX, INITIAL_CONDITIONS
from julesf.triffid.equations import (
    biomass_pools, litterfall_rate, partitioning_function, triffid_rhs
)

def diagnose_triffid(t_span=(0,50), n_pts=501):
    # 1) run the model
    params = {**PFT_PARAMS, 'competition': COMPETITION_MATRIX}
    t, Lb, nu = solve_triffid(t_span, INITIAL_CONDITIONS, params, method='RK45')

    # 2) make Lb “safe” for any negative‐exponent: floor at tiny positive
    Lb_safe = np.maximum(Lb, 1e-8)

    # 3) reshape all PFT parameter vectors to (n_pfts,1)
    def vec(key): 
        return np.array([PFT_PARAMS[p][key] for p in PFT_PARAMS])[:, None]

    sigma1 = vec('sigma1')
    awl    = vec('awl')
    bwl    = vec('bwl')
    gamma_l = vec('gamma_l')
    gamma_r = vec('gamma_r')
    gamma_w = vec('gamma_w')
    Lmin   = vec('Lmin')
    Lmax   = vec('Lmax')

    # 4) pools & fluxes
    L, R, W_pool, Cv = biomass_pools(Lb_safe, {'sigma1': sigma1,
                                               'awl':    awl,
                                               'bwl':    bwl})
    lam    = partitioning_function(Lb_safe, Lmin, Lmax)
    Lambda1 = litterfall_rate(L, R, W_pool, {'gamma_l': gamma_l,
                                             'gamma_r': gamma_r,
                                             'gamma_w': gamma_w})

    # 5) compute rates of change from the TRIFFID RHS
    n_pfts = Lb_safe.shape[0]
    # stack dy/dt for each timestep
    dy = np.vstack([
        triffid_rhs(
            t_i,
            np.concatenate([Lb[:, idx], nu[:, idx]]),
            params
        )
        for idx, t_i in enumerate(t)
    ])  # shape (nt, 2*n_pfts)

    # reshape into (n_pfts, nt)
    dLb_dt = dy[:, :n_pfts].T
    dnu_dt = dy[:, n_pfts:].T

    return dict(
        t=t,
        Lb=Lb,
        nu=nu,
        L=L,
        R=R,
        W=W_pool,
        Cv=Cv,
        lam=lam,
        Lambda1=Lambda1,
        dLb_dt=dLb_dt,
        dnu_dt=dnu_dt
    )

def plot_triffid_diagnostics(res):
    t     = res['t']
    Lb    = res['Lb']
    nu    = res['nu']
    L, R, W, Cv = res['L'], res['R'], res['W'], res['Cv']
    lam   = res['lam']
    Lambda1 = res['Lambda1']
    dLb_dt  = res['dLb_dt']
    dnu_dt  = res['dnu_dt']
    pfts = list(PFT_PARAMS.keys())

    fig, axs = plt.subplots(4, 2, figsize=(12, 10), sharex=True)
    # row1: ν
    for i,ax in enumerate(axs[0]):
        ax.plot(t, nu[i], label=f'ν {pfts[i]}')
    axs[0,0].plot(t, nu.sum(axis=0), 'k--', label='∑ν')
    axs[0,0].legend(); axs[0,0].set_ylabel('ν')
    # row2: Lb and dLb/dt
    for i,ax in enumerate(axs[1]):
        ax.plot(t, Lb[i], label=f'Lᵦ {pfts[i]}')
        ax.plot(t, dLb_dt[i], '--', label='dLᵦ/dt')
        ax.legend(); ax.set_ylabel('Lᵦ, dLᵦ/dt')
    # row3: pools L,R,W
    for i,ax in enumerate(axs[2]):
        ax.plot(t, L[i], label='L'); ax.plot(t, R[i], label='R'); ax.plot(t, W[i], label='W')
        ax.legend(); ax.set_ylabel('kg C m⁻²')
    # row4: Cv, λ, Λ₁
    for i,ax in enumerate(axs[3]):
        ax.plot(t, Cv[i],    label='Cᵥ')
        ax.plot(t, lam[i],   label='λ')
        ax.plot(t, Lambda1[i], label='Λ₁')
        ax.legend(); ax.set_ylabel('Cᵥ / λ / Λ₁')

    for ax in axs.flatten():
        ax.grid(alpha=0.2)
    axs[-1,0].set_xlabel('Time (yr)')
    plt.tight_layout()
    plt.show()
    return fig

if __name__ == '__main__':
    res = diagnose_triffid()
    plot_triffid_diagnostics(res)