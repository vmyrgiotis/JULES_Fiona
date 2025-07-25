import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path

from ebm.run_ebm import run_ebm, drivers as ebm_drivers
from physiology.simulation import generate_forcings, run_npp
from physiology.parameters import SIM_SETTINGS


def main():
    # Time settings: 1-week at 0.5h timestep
    t_hours = np.arange(0, 7 * 24 + 0.5, 0.5)       # 0 to 168 h
    dt_days = 0.5 / 24.0                            # 0.5 h in days
    t_span = (0.0, 7.0)                             # days

    # Run EBM
    T0 = ebm_drivers['Tair'](0)                     # initial surface temp in K
    t_ebm, Ts_ebm = run_ebm(t_span, T0, ebm_drivers, dt_out=dt_days)
    t_h = t_ebm * 24.0                              # convert days to hours

    # Prepare physiology forcings
    # Override generate_forcings so T array comes from EBM
    def gen_forcings_override(settings):
        # settings used for dt_hours, days, but ignore built-in T
        days = settings['days']
        dt = settings['dt_hours']
        # build time array
        t = np.arange(0, days * 24, dt)
        # T in °C for physiology: convert EBM Kelvin to Celsius
        T_leaf_C = Ts_ebm - 273.15
        # reuse other drivers from original generate_forcings for I_par, ci, O2
        _, _, I_par, ci, O2 = generate_forcings(settings)
        return t, T_leaf_C, I_par, ci, O2

    import physiology.simulation as sim
    sim.generate_forcings = gen_forcings_override

    # Run physiology NPP with overridden temperature
    t_npp, npp = run_npp('needleleaf_tree')
    # t_npp corresponds to hours since 0 with dt_hours = 0.5

    # Plot results
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    ax1.plot(t_h, Ts_ebm, label='Surface T* (K)')
    ax1.set_ylabel('Temperature (K)')
    ax1.legend()

    ax2.plot(t_npp, npp['Pi_net'], label='NPP')
    ax2.set_ylabel('Net Primary Productivity')
    ax2.set_xlabel('Time (hours)')
    ax2.legend()

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    # Update SIM_SETTINGS for 0.5h, 7 days
    SIM_SETTINGS['days'] = 7
    SIM_SETTINGS['dt_hours'] = 0.5
    main()
