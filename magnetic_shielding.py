# Magnetic shielding for the case of a cylindrical shell and for the case of a spherical shell.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
import math

# =================================================================
# ParametersInput parameters
# =================================================================

mu_r = 1000.0  # Relative permeability (typical for ferromagnetic materials)
H0 = 1.0       # Intensity of the applied magnetic field
a = 1.0        # Inner Radius (m)
b = 2.0        # External Radius (m)
d = b - a      # Layer thickness
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
    
    # Field in the internal cavity
    H_int = (9 * mu_r * H0) / denominator
    
    # Shielding Factor
    SF = H0 / H_int
    
    # Constants for plotting
    A = H0 * (3 * (2 * mu_r + 1)) / denominator
    B = H0 * (3 * (mu_r - 1) * a**3) / denominator
    C = H0 * b**3 * ((mu_r - 1) * (2 * mu_r + 1) - 2 * (a**3 / b**3) * (mu_r - 1)**2) / denominator
    
    return H_int, SF, A, B, C

H_int_esf, SF_esf, A_esf, B_esf, C_esf = calculate_sphere_shielding(H0, mu_r, a, b)

n = 100
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

# 1. Cavity region (r < a)
Hx_esf[cavity] = H_int_esf
Hy_esf[cavity] = 0.0

# 2. Region with magnetic material (a ≤ r ≤ b)
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
    
    # Field in the internal cavity (H_int)
    H_int = (4 * mu_r * b**2) / denom * H0
    
    # Shielding Factor
    SF = H0 / H_int
    
    # Constantes para plotagem [8]
    A = 2 * (mu_r + 1) * b**2 / denom
    B = 2 * (mu_r - 1) * a**2 * b**2 / denom
    C = b**2 * (mu_r**2 - 1) * (b**2 - a**2) / denom
    
    return H_int, SF, A, B, C

H_int_cil, SF_cil, A_cil, B_cil, C_cil = calculate_cylinder_shielding(H0, mu_r, a, b)


Hx_cil = np.zeros_like(X)
Hy_cil = np.zeros_like(Y)

# 1. Cavity region (r < a)
Hx_cil[cavity] = H_int_cil
Hy_cil[cavity] = 0.0

# 2. Layer Region (a ≤ r ≤ b)
r_layer = R[layer]
theta_layer = Theta[layer]
Hr_layer = (A_cil - B_cil / r_layer**2) * H0 * np.cos(theta_layer)
Htheta_layer = -(A_cil + B_cil / r_layer**2) * H0 * np.sin(theta_layer)
Hx_cil[layer] = Hr_layer * np.cos(theta_layer) - Htheta_layer * np.sin(theta_layer)
Hy_cil[layer] = Hr_layer * np.sin(theta_layer) + Htheta_layer * np.cos(theta_layer)

# 3. External region (r > b)
r_ex = R[exterior]
theta_ex = Theta[exterior]
Hr_ex = H0 * (1 + C_cil / r_ex**2) * np.cos(theta_ex)
Htheta_ex = -H0 * (1 - C_cil / r_ex**2) * np.sin(theta_ex)
Hx_cil[exterior] = Hr_ex * np.cos(theta_ex) - Htheta_ex * np.sin(theta_ex)
Hy_cil[exterior] = Hr_ex * np.sin(theta_ex) + Htheta_ex * np.cos(theta_ex)


# =================================================================
# Results View
# =================================================================

fig, axes = plt.subplots(1, 2, figsize=(18, 8), dpi=100)

# Plotting the SPHERE
ax1 = axes[0]
ax1.set_title(f'(a) Magnetic Shielding - Sphere ($\\mu_r$ = {mu_r})', fontsize=20, fontweight='bold')
ax1.streamplot(X, Y, Hx_esf, Hy_esf, density=[1.0, 2.0], color='blue', linewidth=1.2, arrowsize=1.7)
ax1.add_patch(Circle((0, 0), b, color='lightgray', alpha=0.3, zorder=2))
ax1.add_patch(Circle((0, 0), a, color='white', ec='k', linewidth=1.5, zorder=3))
ax1.set_xlim(-L, L)
ax1.set_ylim(-L, L)
ax1.set_xlabel('X (m)', fontsize=20, fontweight='bold')
ax1.set_ylabel('Y (m)', fontsize=20, fontweight='bold')
ax1.set_aspect('equal', adjustable='box')
#ax1.text(0.1, 0.5, f'$H_{{int}} = {H_int_esf:.4f}$', transform=ax1.transAxes, fontsize=12, bbox=dict(facecolor='white', alpha=0.8))


# Cylinder plot
ax2 = axes[1]
ax2.set_title(f'(b) Magnetic Shielding - Cylinder ($\\mu_r$ = {mu_r})', fontsize=20, fontweight='bold')
ax2.streamplot(X, Y, Hx_cil, Hy_cil, density=[1.0, 2.0], color='green', linewidth=1.2, arrowsize=1.7)
ax2.add_patch(Circle((0, 0), b, color='lightgray', alpha=0.3, zorder=2))
ax2.add_patch(Circle((0, 0), a, color='white', ec='k', linewidth=1.5, zorder=3))
ax2.set_xlim(-L, L)
ax2.set_ylim(-L, L)
ax2.set_xlabel('X (m)', fontsize=20, fontweight='bold')
ax2.set_aspect('equal', adjustable='box')
#ax2.text(0.1, 0.5, f'$H_{{int}} = {H_int_cil:.4f}$', transform=ax2.transAxes, fontsize=12, bbox=dict(facecolor='white', alpha=0.8))

plt.tight_layout()
plt.show()

# =================================================================
# Shielding Factor (SF)
# =================================================================

#print("\n" + "="*50)
#print(f"RESULTS FOR SHIELDING (H0 = {H0:.1f}, μr = {mu_r:.1f}, a={a:.1f}, b={b:.1f})")
#print("="*50)

# SPHERE
#print(f"1. SPHERE (Based on the solution for spherical shielding):")
#print(f"   Residual Internal Field (H_int): {H_int_esf:.6f} * H0")
#print(f"   Shielding Factor (SF): {SF_esf:.2f}")

# Cylinder
#print(f"\n2. INFINITE CYLINDER (Field perpendicular to the axis):")
#print(f"   Residual Internal Field (H_int): {H_int_cil:.6f} * H0")
#print(f"   Shielding Factor (SF_T): {SF_cil:.2f}")

#if SF_esf > SF_cil:
#    print(f"In this example, SF_sphere ({SF_esf:.2f}) is larger than SF_cylinder ({SF_cil:.2f}) by a factor of {SF_esf/SF_cil:.2f}.")
