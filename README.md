# Magnetic Shielding Analysis Tool

## Overview

This software package provides computational tools for analyzing magnetic shielding effectiveness in two common geometries: spherical shells and cylindrical shells under transverse magnetic fields. The implementation is based on analytical solutions of magnetostatic boundary value problems, following classical electromagnetic theory.

## System Requirements

### Software Prerequisites
- Python 3.8 or higher
- NumPy (≥ 1.20)
- Matplotlib (≥ 3.5)

### Hardware Requirements
- Minimum 2GB RAM
- 100MB available disk space
- Display capable of 1280×720 resolution

## Installation

### Method 1: Direct Python Environment
# Clone or download the source files
git clone <https://github.com/xedwfe/magneticshielding>
cd magnetic-shielding

# Install required packages
pip install numpy matplotlib

### Method 2: Conda Environment
conda create -n magnetic_shielding python=3.9
conda activate magnetic_shielding
conda install numpy matplotlib

## Files Included

1. magnetic_shielding.py - Primary simulation module for field visualization
2. magnetic_shielding_S_graf.py - Shielding factor analysis module
3. README.md - This documentation file

## Usage Instructions

### 1. Magnetic Field Visualization

To generate magnetic field streamlines for both spherical and cylindrical geometries:

python magnetic_shielding.py

#### Input Parameters (Modifiable in Code)
mu_r = 1000.0    # Relative permeability (dimensionless)
H0 = 1.0         # Applied field intensity (A/m)
a = 1.0          # Inner radius (m)
b = 2.0          # Outer radius (m)
L = 4.0          # Plot domain limit (m)

#### Output
- Two-panel figure showing magnetic field streamlines
- Left panel: Spherical shell shielding
- Right panel: Cylindrical shell shielding
- Visual representation of field attenuation in cavity regions

### 2. Shielding Factor Analysis

To analyze shielding factor dependence on material permeability:

python magnetic_shielding_S_graf.py

#### Customizable Parameters
mu_r_min = 1           # Minimum relative permeability (log scale start)
mu_r_max = 50000       # Maximum relative permeability (log scale end)
num_points = 500       # Resolution of permeability sweep

#### Output
- Single log-log plot comparing shielding factors
- Analytical asymptotic approximations
- Performance comparison at key permeability values
- Technical annotations for material classes

## Theory and Equations

### Spherical Shell
Internal Field:
H_int = (9μ_r·H₀) / [(μ_r+2)(2μ_r+1) - 2(a³/b³)(μ_r-1)²]

Shielding Factor (asymptotic):
SF_sphere ≈ (2/9)μ_r[1 - (a/b)³] + 1 (for μ_r ≫ 1)

### Cylindrical Shell (Transverse Field)
Internal Field:
H_int = (4μ_r·b²·H₀) / [(μ_r+1)²b² - (μ_r-1)²a²]

Shielding Factor (asymptotic):
SF_cylinder ≈ (1/4)μ_r[1 - (a/b)²] + 1 (for μ_r ≫ 1)

## Parameter Guidelines

### Material Properties
Material Class     Typical μ_r Range   Application
Air/Vacuum        1.0                 Reference
Ferrites          100-1000            Low-frequency shielding
Silicon Steel     1000-5000           Power transformers
Ni-Fe Alloys      5000-100000         Precision shielding

### Geometric Considerations
- Aspect ratio (b/a) typically ranges 1.1-2.0
- Thickness d = b-a affects shielding effectiveness
- Larger d/(2b) ratios improve shielding but increase mass

## Interpretation of Results

### Field Plots
- Streamline density indicates field intensity
- Arrow direction shows field orientation
- White cavity region displays field attenuation
- Gray annulus represents magnetic material

### Shielding Factor Plot
- Y-axis: Shielding Factor SF = H₀/H_int
- Higher SF values indicate better shielding
- Slope indicates sensitivity to permeability changes
- Performance ratio quantifies geometric advantage

## Troubleshooting

### Common Issues

1. No graphical output appears
   - Verify Matplotlib installation: python -c "import matplotlib.pyplot"
   - Check display backend settings
   - Ensure script execution completes without errors

2. Memory errors with high-resolution meshes
   - Reduce n parameter in magnetic_shielding.py (line 39)
   - Decrease num_points in magnetic_shielding_S_graf.py

3. Unphysical results
   - Verify μ_r > 1 for magnetic materials
   - Ensure b > a > 0
   - Check for division by zero in asymptotic regimes

### Numerical Stability
- Solutions valid for μ_r ≥ 1
- High-μ_r approximations accurate for μ_r > 100
- Geometric singularities avoided by a/b < 1

## Validation

The code implements analytical solutions verified against:
- Jackson, Classical Electrodynamics (3rd ed.), Section 5.12
- Hoburg, J.F. (1995), IEEE Transactions on Magnetics
- Smythe, W.R. (1950), Static and Dynamic Electricity

## Extensions and Modifications

### Adding New Geometries
1. Implement new shielding function following existing pattern
2. Add to plot_shielding_factor_analysis() function
3. Include in comparative visualization

### Custom Output Formats
Modify Matplotlib commands to:
- Change color schemes
- Adjust streamline densities
- Export data for external processing

## License and Attribution

This software is provided for educational and research purposes. Users are encouraged to cite relevant electromagnetic theory references when publishing results obtained with this tool.

## Support

For technical questions or bug reports, please provide:
1. Complete error messages
2. Parameter settings used
3. Python environment details

---
Version: 1.0 | Last Updated: December 2025
