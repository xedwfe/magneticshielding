# Magnetic Shielding Analysis Tool

<!-- [CONFIRMAR]: ajuste "xedwfe/magneticshielding" e o branch "main" nos dois
     badges abaixo para o nome exato do repositório antes do commit. -->
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/xedwfe/magneticshielding/blob/main/magnetic_shielding_interactive.ipynb)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/xedwfe/magneticshielding/main?filepath=magnetic_shielding_interactive.ipynb)
![Python](https://img.shields.io/badge/python-3.8%2B-blue)

## Overview

This software package provides computational tools for analyzing magnetic
shielding effectiveness in two common geometries: spherical shells and
cylindrical shells under transverse magnetic fields. The implementation is
based on **exact analytical solutions** of magnetostatic boundary value
problems, following classical electromagnetic theory.

It is the companion code of the paper *"Magnetic shielding as a
motivational tool in teaching classical electromagnetic theory"*
(submitted to the *European Journal of Physics*, 2026).
<!-- [CONFIRMAR]: adicionar DOI do artigo e do release Zenodo após o aceite. -->

## Quick start — no installation required

The interactive notebook runs in any web browser, with sliders for
μᵣ, a/b and the remaining parameters:

- **Google Colab:** click the *Open in Colab* badge above (requires only a
  Google account).
- **Binder:** click the *Binder* badge above (no account at all; first
  launch may take a minute to build).

This zero-installation route is intended for classrooms and teaching
laboratories where installing software is restricted.

## Local installation

### Software prerequisites
- Python 3.8 or higher
- NumPy (≥ 1.20)
- Matplotlib (≥ 3.5)
- ipywidgets (≥ 8.0) — only for the interactive notebook

```bash
# Clone or download the source files
git clone https://github.com/xedwfe/magneticshielding.git
cd magneticshielding

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
| `magnetic_shieldingV3_corrected.py` | Field-line visualization (sphere + transverse cylinder), Figure 1 of the paper |
| `magnetic_shielding_S_graf_corrected.py` | Exact shielding-factor analysis, Figure 2 of the paper |
| `magnetic_shielding_interactive.ipynb` | Interactive notebook (sliders; Colab/Binder-ready) |
| `verify_cavity_uniformity.py` | Numerical verification that the cavity field is exactly uniform |
| `requirements.txt` | Dependencies for pip / Binder |
| `README.md` | This documentation file |

<!-- [CONFIRMAR]: se preferir, renomeie os dois scripts retirando o sufixo
     "_corrected" (e atualize a tabela acima e os comandos abaixo). -->

## Usage instructions

### 1. Magnetic field visualization

```bash
python magnetic_shieldingV3_corrected.py
```

**Input parameters (modifiable in code)**

```python
mu_r = 1000.0    # Relative permeability (dimensionless)
H0 = 1.0         # Applied field intensity (A/m)
a = 1.0          # Inner radius (m)
b = 2.0          # Outer radius (m)
L = 4.0          # Plot domain limit (m)
```

**Output** — a two-panel figure (`Figure_1_modelstyle.png`): spherical
shell (left, meridian plane) and cylindrical shell (right, transverse
plane), with the streamlines concentrated in the material and three
straight lines showing the small, exactly uniform residual field inside
each cavity. At a line density proportional to the flux, no line would
cross the cavity at μᵣ = 1000 — the interior line count is illustrative.

### 2. Shielding factor analysis

```bash
python magnetic_shielding_S_graf_corrected.py
```

**Customizable parameters**

```python
mu_r_min = 1           # Minimum relative permeability (log scale start)
mu_r_max = 50000       # Maximum relative permeability (log scale end)
num_points = 500       # Resolution of permeability sweep
```

**Output** — a log–log plot (`Figure_2_modelstyle.png`) of the **exact**
shielding factor for both geometries, the high-permeability asymptotes
(dashed), the no-shielding reference SF = 1, and physically correct
material annotations.

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

$$H_{int} = \frac{9\,\mu_r H_0}{(\mu_r+2)(2\mu_r+1) - 2(a^3/b^3)(\mu_r-1)^2}$$

Shielding factor — **exact closed form**, with $k_3 = (a/b)^3$:

$$SF_{sph} = \frac{2}{9}(1-k_3)\,\mu_r + \frac{5+4k_3}{9} + \frac{2(1-k_3)}{9\mu_r}$$

### Cylindrical shell (transverse field)

Internal field:

$$H_{int} = \frac{4\,\mu_r b^2 H_0}{(\mu_r+1)^2 b^2 - (\mu_r-1)^2 a^2}$$

Shielding factor — **exact closed form**, with $k_2 = (a/b)^2$:

$$SF_{cyl} = \frac{1}{4}(1-k_2)\,\mu_r + \frac{1+k_2}{2} + \frac{1-k_2}{4\mu_r}$$

Both exact expressions reduce to SF = 1 at μᵣ = 1 (no material, no
shielding). For μᵣ ≫ 1 they approach the asymptotes
$\frac{2}{9}\mu_r(1-k_3)$ and $\frac{1}{4}\mu_r(1-k_2)$ — note that the
commonly quoted "+1" in these asymptotic formulas is an ad hoc addition
that overestimates SF by 19–27% for μᵣ ≤ 10; the exact forms above should
be preferred. The crossover between the low-effectiveness plateau
(SF ≈ 1) and the linear-growth regime occurs at
μᵣ\* ~ 9/[2(1−k₃)] (sphere) and 4/(1−k₂) (cylinder).

### Why is the cavity field uniform?

Within the ideal model the interior field is *exactly* uniform: the
interior solution regular at the origin can only contain terms
rˡPₗ(cos θ), a uniform applied field excites only ℓ = 1, and every other
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
  dramatic geometric contrast is transverse vs **axial**: an infinite open
  cylinder provides no static shielding along its axis.

## Interpretation of results

### Field plots
- Streamline density outside/inside the material indicates field intensity.
- Arrow direction shows field orientation.
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
  Compatibility* **37**(4), 574–579, doi:10.1109/15.477342;
- Smythe, W. R. (1950), *Static and Dynamic Electricity*.

All boundary conditions are satisfied to machine precision (residuals
~10⁻¹⁶), and the interior-field uniformity is verified independently by
`verify_cavity_uniformity.py`.

## Changelog

### v2.0 (July 2026)
- **Bug fix (sphere, exterior region):** the induced-dipole coefficient
  was transcribed incorrectly; the corrected coefficient restores the
  boundary conditions at r = b to machine precision.
- Shielding-factor plot now uses the **exact closed forms** (the previous
  version plotted the "+1" high-permeability approximation as if exact).
- Field lines inside the cavity are now visible (previously hidden by the
  cavity patch), drawn as three straight seeded lines.
- Material annotations corrected (air/polymers have μᵣ ≈ 1, not 10).
- Hoburg (1995) reference corrected (IEEE Trans. Electromagn. Compat.,
  not IEEE Trans. Magn.).
- Added: interactive Colab/Binder notebook and numerical verification
  script.

### v1.0 (December 2025)
- Initial release.

## How to cite

If you use this code, please cite the companion paper:

> M. S. Marques, R. R. Marques, J. J. da Silva and E. F. de Almeida
> Júnior, "Magnetic shielding as a motivational tool in teaching
> classical electromagnetic theory", submitted to *Eur. J. Phys.* (2026).
<!-- [CONFIRMAR]: atualizar com volume/páginas/DOI após publicação, e com o
     DOI Zenodo do release arquivado do código. -->

## Transparency note

An AI assistant (Claude, Anthropic) was used to help draft and refactor
parts of this code and documentation; all physics, derivations and
results were reviewed, tested and verified by the authors, who take full
responsibility for the content.

## License and attribution

This software is provided for educational and research purposes.
<!-- [CONFIRMAR]: recomendo adicionar um arquivo LICENSE explícito (MIT é a
     escolha usual para código educacional) e ajustar esta seção. -->

## Support

For technical questions or bug reports, please open an issue and provide:
1. Complete error messages
2. Parameter settings used
3. Python environment details

---
Version: 2.0 | Last updated: July 2026
# Magnetic shielding for the case of a cylindrical shell and for the case
# of a spherical shell.
#
# CORRECTED VERSION of magnetic_shieldingV3.py -- same visual style as the
# original model. Changes:
#
#   (1) BUG FIX (sphere, exterior region): the induced-dipole coefficient
#       C_esf was transcribed incorrectly from the Supplementary Material.
#       Correct numerator (coefficient G1 of the SM):
#           (mu_r - 1)(2 mu_r + 1)(1 - a^3/b^3)
#       The wrong C violated the boundary conditions at r = b (relative
#       error ~9% at mu_r = 2). With the fix, continuity of B_n and H_t
#       holds to machine precision (see verify.py).
#   (2) Field lines are now VISIBLE inside the cavity (Reviewer 2,
#       suggestion 1). The white cavity patch no longer covers the lines
#       (zorder), the cavity is excluded from the density-based streamplot
#       (which would fill it with as many lines as the exterior, visually
#       overstating the weak interior field), and exactly THREE seeded
#       lines are drawn inside.
#       PHYSICS NOTE: within this model the cavity field is EXACTLY
#       uniform, H_vec = H_int * x_hat. The lines inside are therefore
#       STRAIGHT, horizontal and equally spaced -- a small residual field
#       does NOT imply curvature. Only the l = 1 (nu = 1) harmonic is
#       excited by a uniform applied field; every other harmonic satisfies
#       a homogeneous system with nonzero determinant and vanishes
#       (see verify_cavity_uniformity.py for the numerical demonstration).
#   (3) H0 harmonized: all coefficients now carry H0 (the cylinder block
#       used dimensionless coefficients multiplied by H0 afterwards).
#   (4) Cleanup: unused imports (Rectangle, math) and unused d removed;
#       stray "[8]" comment marker removed; finer grid (n = 300) for
#       smoother streamlines.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# =================================================================
# Input parameters
# =================================================================

mu_r = 1000.0  # Relative permeability (typical for ferromagnetic materials)
H0 = 1.0       # Intensity of the applied magnetic field
a = 1.0        # Inner Radius (m)
b = 2.0        # External Radius (m)
L = 4.0        # Plot limit

# -----------------------------------------------------------------
# 1. CALCULATION AND SIMULATION FOR SPHERES
# -----------------------------------------------------------------

def calculate_sphere_shielding(H0, mu_r, a, b):
    """
    Calculate the internal magnetic field (H_int) and the Shielding Factor (SF)
    for a spherical shell under a uniform field H0.
    """

    # General Denominator
    denominator = (mu_r + 2) * (2 * mu_r + 1) - 2 * (a**3 / b**3) * (mu_r - 1)**2

    # Field in the internal cavity (exactly uniform: H = H_int * x_hat)
    H_int = (9 * mu_r * H0) / denominator

    # Shielding Factor
    SF = H0 / H_int

    # Constants for plotting
    A = H0 * (3 * (2 * mu_r + 1)) / denominator
    B = H0 * (3 * (mu_r - 1) * a**3) / denominator
    # --- corrected exterior (induced dipole) coefficient: C = G1 of the SM ---
    C = H0 * b**3 * (mu_r - 1) * (2 * mu_r + 1) * (1 - a**3 / b**3) / denominator

    return H_int, SF, A, B, C

H_int_esf, SF_esf, A_esf, B_esf, C_esf = calculate_sphere_shielding(H0, mu_r, a, b)

n = 300
x = np.linspace(-L, L, n)
y = np.linspace(-L, L, n)
X, Y = np.meshgrid(x, y)
R = np.sqrt(X**2 + Y**2)
Theta = np.arctan2(Y, X)

Hx_esf = np.zeros_like(X)
Hy_esf = np.zeros_like(Y)

cavity = R < a
layer = (R >= a) & (R <= b)
exterior = R > b

# 1. Cavity region (r < a): exactly uniform field along +x
Hx_esf[cavity] = H_int_esf
Hy_esf[cavity] = 0.0

# 2. Region with magnetic material (a <= r <= b)
r_layer = R[layer]
theta_layer = Theta[layer]
Hr_layer = (A_esf - 2 * B_esf / r_layer**3) * np.cos(theta_layer)
Htheta_layer = -(A_esf + B_esf / r_layer**3) * np.sin(theta_layer)
Hx_esf[layer] = Hr_layer * np.cos(theta_layer) - Htheta_layer * np.sin(theta_layer)
Hy_esf[layer] = Hr_layer * np.sin(theta_layer) + Htheta_layer * np.cos(theta_layer)

# 3. External region (r > b)
r_ex = R[exterior]
theta_ex = Theta[exterior]
Hr_ex = (H0 + 2 * C_esf / r_ex**3) * np.cos(theta_ex)
Htheta_ex = -(H0 - C_esf / r_ex**3) * np.sin(theta_ex)
Hx_esf[exterior] = Hr_ex * np.cos(theta_ex) - Htheta_ex * np.sin(theta_ex)
Hy_esf[exterior] = Hr_ex * np.sin(theta_ex) + Htheta_ex * np.cos(theta_ex)


# -----------------------------------------------------------------
# 2. CALCULATION AND SIMULATION FOR CYLINDER (TRANSVERSE FIELD)
# -----------------------------------------------------------------

def calculate_cylinder_shielding(H0, mu_r, a, b):
    """
    Calculate the internal magnetic field (H_int) and the Shielding Factor (SF)
    for a cylindrical shell under a static transverse field H0.
    """

    # Denominator of the analytical formula (Transverse Cylinder)
    denom = (mu_r + 1)**2 * b**2 - (mu_r - 1)**2 * a**2

    # Field in the internal cavity (exactly uniform: H = H_int * x_hat)
    H_int = (4 * mu_r * b**2) / denom * H0

    # Shielding Factor
    SF = H0 / H_int

    # Constants for plotting (H0 included, as in the sphere block)
    A = 2 * (mu_r + 1) * b**2 * H0 / denom
    B = 2 * (mu_r - 1) * a**2 * b**2 * H0 / denom
    C = b**2 * (mu_r**2 - 1) * (b**2 - a**2) * H0 / denom

    return H_int, SF, A, B, C

H_int_cil, SF_cil, A_cil, B_cil, C_cil = calculate_cylinder_shielding(H0, mu_r, a, b)


Hx_cil = np.zeros_like(X)
Hy_cil = np.zeros_like(Y)

# 1. Cavity region (r < a): exactly uniform field along +x
Hx_cil[cavity] = H_int_cil
Hy_cil[cavity] = 0.0

# 2. Layer Region (a <= r <= b)
r_layer = R[layer]
theta_layer = Theta[layer]
Hr_layer = (A_cil - B_cil / r_layer**2) * np.cos(theta_layer)
Htheta_layer = -(A_cil + B_cil / r_layer**2) * np.sin(theta_layer)
Hx_cil[layer] = Hr_layer * np.cos(theta_layer) - Htheta_layer * np.sin(theta_layer)
Hy_cil[layer] = Hr_layer * np.sin(theta_layer) + Htheta_layer * np.cos(theta_layer)

# 3. External region (r > b)
r_ex = R[exterior]
theta_ex = Theta[exterior]
Hr_ex = (H0 + C_cil / r_ex**2) * np.cos(theta_ex)
Htheta_ex = -(H0 - C_cil / r_ex**2) * np.sin(theta_ex)
Hx_cil[exterior] = Hr_ex * np.cos(theta_ex) - Htheta_ex * np.sin(theta_ex)
Hy_cil[exterior] = Hr_ex * np.sin(theta_ex) + Htheta_ex * np.cos(theta_ex)


# =================================================================
# Results View
# =================================================================

def plot_panel(ax, Hx, Hy, color, title):
    """One panel in the style of the original model, with the cavity lines
    visible: the density-based streamplot is applied OUTSIDE the cavity
    only, and three seeded straight lines represent the (exactly uniform)
    interior field."""
    # patches BELOW the lines, so that the cavity no longer hides them
    ax.add_patch(Circle((0, 0), b, color='lightgray', alpha=0.3, zorder=1))
    ax.add_patch(Circle((0, 0), a, facecolor='white', ec='k',
                        linewidth=1.5, zorder=1))

    # density-based lines everywhere EXCEPT the cavity
    Hx_out = np.ma.array(Hx, mask=cavity)
    Hy_out = np.ma.array(Hy, mask=cavity)
    ax.streamplot(X, Y, Hx_out, Hy_out, density=[1.0, 2.0], color=color,
                  linewidth=1.2, arrowsize=1.7)

    # three seeded lines INSIDE the cavity. The interior field is exactly
    # uniform (see header note), so these are straight horizontal lines;
    # their number is illustrative -- at true flux scaling (line density
    # proportional to B) no line would cross the cavity at mu_r = 1000.
    Hx_in = np.ma.array(Hx, mask=~cavity)
    Hy_in = np.ma.array(Hy, mask=~cavity)
    seeds = np.column_stack([np.zeros(3), np.linspace(-0.55 * a, 0.55 * a, 3)])
    try:
        ax.streamplot(X, Y, Hx_in, Hy_in, start_points=seeds, color=color,
                      linewidth=1.0, arrowsize=1.4, broken_streamlines=False)
    except TypeError:  # matplotlib < 3.6
        ax.streamplot(X, Y, Hx_in, Hy_in, start_points=seeds, color=color,
                      linewidth=1.0, arrowsize=1.4)

    ax.set_xlim(-L, L)
    ax.set_ylim(-L, L)
    ax.set_xlabel('X (m)', fontsize=20, fontweight='bold')
    ax.set_aspect('equal', adjustable='box')
    ax.set_title(title, fontsize=20, fontweight='bold')


fig, axes = plt.subplots(1, 2, figsize=(18, 8), dpi=100)

# Plotting the SPHERE
plot_panel(axes[0], Hx_esf, Hy_esf, 'blue',
           f'(a) Magnetic Shielding - Sphere ($\\mu_r$ = {mu_r:g})')
axes[0].set_ylabel('Y (m)', fontsize=20, fontweight='bold')

# Cylinder plot
plot_panel(axes[1], Hx_cil, Hy_cil, 'green',
           f'(b) Magnetic Shielding - Cylinder ($\\mu_r$ = {mu_r:g})')

plt.tight_layout()
plt.savefig('Figure_1_modelstyle.png', bbox_inches='tight')
plt.show()

# =================================================================
# Shielding Factor (SF)
# =================================================================

print("\n" + "=" * 50)
print(f"RESULTS (H0 = {H0:.1f}, mu_r = {mu_r:.1f}, a = {a:.1f}, b = {b:.1f})")
print("=" * 50)
print(f"1. SPHERE   : H_int = {H_int_esf:.6f} * H0   SF = {SF_esf:.2f}")
print(f"2. CYLINDER : H_int = {H_int_cil:.6f} * H0   SF = {SF_cil:.2f}")
print(f"   SF_sphere / SF_cylinder = {SF_esf / SF_cil:.3f}")
