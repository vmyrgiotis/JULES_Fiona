# run_triffid.py
from julesf.triffid.parameters import PFT_PARAMS, COMPETITION_MATRIX, INITIAL_CONDITIONS
from julesf.triffid.simulation import solve_triffid

import matplotlib.pyplot as plt

def main():
    # Combine parameters
    params = {**PFT_PARAMS, 'competition': COMPETITION_MATRIX}
    
    # Solve
    t_span = (0, 50)
    t, Lb, nu = solve_triffid(t_span, INITIAL_CONDITIONS, params)
    
    # Plot (same as current visualization)
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    axes[0].plot(t, nu[0], label='Tree cover')
    axes[0].plot(t, nu[1], label='Grass cover')
    axes[0].set_ylabel('Fractional cover')
    axes[0].legend()

    axes[1].plot(t, Lb[0], label='Tree Lb')
    axes[1].plot(t, Lb[1], label='Grass Lb')
    axes[1].set_ylabel('Balanced LAI')
    axes[1].set_xlabel('Time (yr)')
    axes[1].legend()

    plt.tight_layout()
    plt.show()
    
    return t, Lb, nu

if __name__ == "__main__":
    results = main()