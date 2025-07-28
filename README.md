# JULES Prototype

This folder contains a lightweight Python implementation of the Joint UK Land Environment Simulator (JULES) developed by Fiona Fang.

# Dependencies

* Python 3.7+
* NumPy
* Matplotlib
* SciPy

# Installation

1. Clone the repository:

   ```bash
   git clone https://fionazfang/JULES_Fiona.git
   ```
2. Create a virtual environment and install dependencies:

   ```bash
   python -m venv venv
   source venv/bin/activate     # Unix/macOS
   venv\Scripts\activate        # Windows
   pip install numpy matplotlib scipy
   ```

# TOPMODEL Prototype

## Structure

```
topmodel/  
├── equation.py        # Core physics routines (transmissivity, baseflow, etc.)  
├── simulation.py      # Simulation driver: sets up inputs, runs time-stepping  
├── visualization.py   # Visualization functions (plots of runoff and intermediates)  
└── run_topmodel.py    # Entry-point script: runs simulation and generates plots
```

## Usage

Run the model and visualize results by executing:

```bash
python run_topmodel.py
```

This script will:

1. Simulate a 5×5 grid over one year (daily time-step) of baseflow (`R_b`) and saturation-excess runoff (`R_se`).
2. Produce two figures:

   * Time series of `R_b` and `R_se` for selected cells.
   * Catchment-mean critical topographic index (`\lambda_c`) and saturated fraction (`f_sat`).

#  Vegetation Coupler - TRIFFID, RothC and physiology

---

*Developed by Fiona Fang.*
