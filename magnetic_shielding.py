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
