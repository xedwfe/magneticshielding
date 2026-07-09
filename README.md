# Magnetic Shielding Analysis Tool

## Overview

This software package provides computational tools for analyzing magnetic
shielding effectiveness in two common geometries: spherical shells and
cylindrical shells under transverse magnetic fields. The implementation is
based on exact analytical solutions of magnetostatic boundary value
problems, following classical electromagnetic theory.

It is the companion code of the paper *"Magnetic shielding as a
motivational tool in teaching classical electromagnetic theory"*
(submitted to the *European Journal of Physics*, 2026).

## Quick start — no installation required

The interactive notebook `magnetic_shielding_interactive.ipynb` runs in a
web browser, with sliders for μᵣ, a/b and the remaining parameters, and
requires no local Python installation:

- **Google Colab:** open Colab, choose *File → Open notebook → GitHub*,
  paste this repository's address and select the notebook (a free Google
  account is enough); or download the notebook file and use
  *File → Upload notebook*.
- **Reading only:** GitHub renders the notebook directly in the repository
  page, including the saved figures.

This zero-installation route is intended for classrooms and teaching
laboratories where installing software is restricted.

## System requirements

- Python 3.8 or higher
- NumPy (≥ 1.20)
- Matplotlib (≥ 3.5)
- ipywidgets (≥ 8.0) — only for the interactive notebook

## Local installation

```bash
# Clone or download the source files
git clone <repository-url>
cd <repository-folder>

# Install required packages
pip install -r requirements.txt
```

Or with conda:

```bash
conda create -n magnetic_shielding python=3.9
conda activate magnetic_shielding
conda install numpy matplotlib ipywidgets
```

## Files included

| File | Purpose |
| --- | --- |
| `magnetic_shieldingV3_corrected.py` | Field-line visualization (sphere + transverse cylinder) — Figure 1 of the paper |
| `magnetic_shielding_S_graf_corrected.py` | Exact shielding-factor analysis — Figure 2 of the paper |
| `magnetic_shielding_interactive.ipynb` | Interactive notebook with sliders |
| `verify_cavity_uniformity.py` | Numerical verification that the cavity field is exactly uniform |
| `requirements.txt` | Dependencies for pip |
| `README.md` | This documentation file |

## Usage instructions

### 1. Magnetic field visualization

```bash
python magnetic_shieldingV3_corrected.py
```

Input parameters (modifiable in code):

```python
mu_r = 1000.0    # Relative permeability (dimensionless)
H0 = 1.0         # Applied field intensity (A/m)
a = 1.0          # Inner radius (m)
b = 2.0          # Outer radius (m)
L = 4.0          # Plot domain limit (m)
```

Output — a two-panel figure (`Figure_1_modelstyle.png`): spherical shell
(left, meridian plane) and cylindrical shell (right, transverse plane),
with the streamlines concentrated in the material and three straight
lines showing the small, exactly uniform residual field inside each
cavity. At a line density proportional to the flux, no line would cross
the cavity at μᵣ = 1000 — the interior line count is illustrative.

### 2. Shielding factor analysis

```bash
python magnetic_shielding_S_graf_corrected.py
```

Customizable parameters:

```python
mu_r_min = 1           # Minimum relative permeability (log scale start)
mu_r_max = 50000       # Maximum relative permeability (log scale end)
num_points = 500       # Resolution of permeability sweep
```

Output — a log–log plot (`Figure_2_modelstyle.png`) of the exact
shielding factor for both geometries, the high-permeability asymptotes
(dashed), the no-shielding reference SF = 1, and material annotations.

### 3. Verification script

```bash
python verify_cavity_uniformity.py
```

Demonstrates, by two independent routes, that the field inside the cavity
is *exactly* uniform (straight field lines): (a) the boundary-condition
systems of all harmonics ℓ ≠ 1 (ν ≠ 1) are homogeneous with non-vanishing
determinants, so only the uniform harmonic is excited; (b) an independent
finite-volume solution of the radial ν = 1 problem gives f(s)/s constant
to machine precision inside the cavity, with clean second-order
convergence of H_int to the exact solution of the same truncated problem.

## Theory and equations

### Spherical shell

Internal field:

    H_int = (9 μr H0) / [(μr + 2)(2μr + 1) − 2(a³/b³)(μr − 1)²]

Shielding factor — exact closed form, with k₃ = (a/b)³:

    SF_sph = (2/9)(1 − k₃) μr + (5 + 4k₃)/9 + 2(1 − k₃)/(9 μr)

### Cylindrical shell (transverse field)

Internal field:

    H_int = (4 μr b² H0) / [(μr + 1)² b² − (μr − 1)² a²]

Shielding factor — exact closed form, with k₂ = (a/b)²:

    SF_cyl = (1/4)(1 − k₂) μr + (1 + k₂)/2 + (1 − k₂)/(4 μr)

Both exact expressions reduce to SF = 1 at μᵣ = 1 (no material, no
shielding). For μᵣ ≫ 1 they approach the asymptotes (2/9)μᵣ(1 − k₃) and
(1/4)μᵣ(1 − k₂) — note that the commonly quoted "+1" in these asymptotic
formulas is an ad hoc addition that overestimates SF by 19–27% for
μᵣ ≤ 10; the exact forms above should be preferred. The crossover between
the low-effectiveness plateau (SF ≈ 1) and the linear-growth regime
occurs at μᵣ* ~ 9/[2(1 − k₃)] (sphere) and 4/(1 − k₂) (cylinder).

### Why is the cavity field uniform?

Within the ideal model the interior field is *exactly* uniform: the
interior solution regular at the origin can only contain terms
rˡ Pₗ(cos θ), a uniform applied field excites only ℓ = 1, and every other
harmonic satisfies a homogeneous system with non-vanishing determinant.
A small residual field is therefore *not* a curved one. Curvature appears
only when the ideal model is broken (applied-field gradients, finite
length or open ends, holes and seams, magnetic saturation). See
`verify_cavity_uniformity.py`.

## Parameter guidelines

### Material properties

| Material class | Typical μᵣ range | Application |
| --- | --- | --- |
| Air / vacuum / polymers | ≈ 1 | No shielding (SF = 1) |
| Ferrites | 10²–10³ | Low-frequency shielding |
| Silicon steel | 10³–5×10³ | Power transformers, typical shielding |
| Ni–Fe alloys (mu-metal) | 5×10³–10⁵ | Precision shielding |

### Geometric considerations

- Aspect ratio (b/a) typically ranges 1.1–2.0.
- In the thin-shell limit, the thickness ratio d/(2b), with d = b − a,
  sets the geometric factor of the shielding; larger ratios improve
  shielding but increase mass.
- For the same a/b, the sphere outperforms the transverse cylinder by at
  most a factor 4/3 (≈ 1.04 at μᵣ = 1000 for a = 1 m, b = 2 m). The
  dramatic geometric contrast is transverse vs axial: an infinite open
  cylinder provides no static shielding along its axis.

## Interpretation of results

### Field plots

- Streamline density outside/inside the material indicates field
  intensity; arrow direction shows field orientation.
- The three straight lines in the white cavity display the direction and
  uniformity of the (small) residual field.
- The gray annulus represents the magnetic material.

### Shielding factor plot

- Y-axis: shielding factor SF = H₀/H_int (exact).
- Higher SF values indicate better shielding.
- The unit slope of the log–log curve beyond the crossover reflects
  SF ∝ μᵣ in the linear-material model; in real materials, saturation
  (μᵣ decreasing with H) sets a practical ceiling.

## Troubleshooting

1. **No graphical output appears** — verify Matplotlib
   (`python -c "import matplotlib.pyplot"`), check the display backend,
   or use the saved PNG files.
2. **Memory errors with high-resolution meshes** — reduce `n` in
   `magnetic_shieldingV3_corrected.py` or `num_points` in
   `magnetic_shielding_S_graf_corrected.py`.
3. **Unphysical results** — ensure μᵣ ≥ 1 and b > a > 0.

## Validation

The code implements analytical solutions verified against:

- Jackson, *Classical Electrodynamics* (3rd ed.), Section 5.12;
- Hoburg, J. F. (1995), *IEEE Transactions on Electromagnetic
  Compatibility* **37**(4), 574–579;
- Smythe, W. R. (1950), *Static and Dynamic Electricity*.

All boundary conditions are satisfied to machine precision (residuals
~10⁻¹⁶), and the interior-field uniformity is verified independently by
`verify_cavity_uniformity.py`.

## Changelog

### v2.0 (July 2026)

- Bug fix (sphere, exterior region): the induced-dipole coefficient was
  transcribed incorrectly; the corrected coefficient restores the
  boundary conditions at r = b to machine precision.
- Shielding-factor plot now uses the exact closed forms (the previous
  version plotted the "+1" high-permeability approximation as if exact).
- Field lines inside the cavity are now visible (previously hidden by the
  cavity patch), drawn as three straight seeded lines.
- Material annotations corrected (air/polymers have μᵣ ≈ 1, not 10).
- Hoburg (1995) reference corrected (IEEE Trans. Electromagn. Compat.,
  not IEEE Trans. Magn.).
- Added: interactive notebook and numerical verification script.

### v1.0 (December 2025)

- Initial release.

## How to cite

If you use this code, please cite the companion paper:

> M. S. Marques, R. R. Marques, J. J. da Silva and E. F. de Almeida
> Júnior, "Magnetic shielding as a motivational tool in teaching
> classical electromagnetic theory", submitted to *Eur. J. Phys.* (2026).

## Transparency note

An AI assistant (Claude, Anthropic) was used to help draft and refactor
parts of this code and documentation; all physics, derivations and
results were reviewed, tested and verified by the authors, who take full
responsibility for the content.

## License and attribution

This software is provided for educational and research purposes. Users
are encouraged to cite the companion paper and the relevant
electromagnetic theory references when publishing results obtained with
this tool.

## Support

For technical questions or bug reports, please open an issue in this
repository and provide:

1. Complete error messages
2. Parameter settings used
3. Python environment details

---
Version: 2.0 | Last updated: July 2026
