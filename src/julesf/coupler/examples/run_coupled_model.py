from coupler.triffid_rothc.coupler import run_coupled_model
from rothc.parameters import POOLS, C0_default
import matplotlib.pyplot as plt

def main():
    # Setup initial conditions
    triffid_init = [6.0, 3.0, 0.8, 0.2]
    rothc_init = [C0_default[p] for p in POOLS]
    
    initial_conditions = {
        'triffid': triffid_init,
        'rothc': rothc_init
    }
    
    drivers_ts = {
        'T_soil': lambda t: 283.15,
        'moisture': lambda t: 0.5,
        's': lambda t: 0.5
    }
    
    # Run model
    results = run_coupled_model(
        t_span=(0, 365),
        initial_conditions=initial_conditions,
        drivers_ts=drivers_ts
    )
    
    # Plot results
    plt.figure(figsize=(10,6))
    plt.plot(results['time'], results['triffid'][2:].T, label=['PFT1', 'PFT2'])
    plt.xlabel('Time (days)')
    plt.ylabel('Vegetation cover')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()