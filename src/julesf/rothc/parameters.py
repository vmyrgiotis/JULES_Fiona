# parameters.py

# Pool names
POOLS = ["DPM", "RPM", "BIO", "HUM"]

# Base respiration rates κ_i [day⁻¹], from table 8
kappa = {
    "DPM": 3.22e-7 * 86400,
    "RPM": 9.65e-9 * 86400,
    "BIO": 2.12e-8 * 86400,
    "HUM": 6.43e-10 * 86400,
}

# Litter partitioning α_dr by PFT
alpha_dr = {
    "tree": 0.25,
    "shrub": 0.33,
    "grass": 0.67,
    "crop": 1.44,
}

# Soil texture and β_R
beta_R = {
    "sand": 0.6,
    "loam": 0.75,
    "clay": 0.85,
    "peat": 0.9,
}
beta_R_default = 0.7

# Default initial conditions
C0_default = {
    "DPM": 0.1,
    "RPM": 0.5,
    "BIO": 0.05,
    "HUM": 1.0,
}
