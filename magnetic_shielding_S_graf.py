"""
Magnetic Shielding Factor Analysis

This code analyzes and plots the shielding factor (SF) for spherical and
cylindrical shells as a function of relative permeability. The analysis is
based on solutions for magnetostatic shielding in hollow geometries.

CORRECTED VERSION of magnetic_shielding_S_graf.py -- same visual style as
the original model. Changes:

  (1) The curves now plot the EXACT shielding factor SF = H0/H_int (the
      exact expression was computed but commented out in the original;
      the plotted quantity was the high-permeability approximation
      "(2/9) mu_r (1 - k3) + 1", whose ad hoc "+1" overestimates SF by
      19-27% for mu_r <= 10 and gives SF(1) = 1.19 instead of 1).
  (2) The exact SF admits a closed form, shown in the equation box:
          SF_sph = (2/9)(1-k3) mu_r + (5+4k3)/9 + 2(1-k3)/(9 mu_r),
          SF_cyl = (1/4)(1-k2) mu_r + (1+k2)/2 + (1-k2)/(4 mu_r),
      with k3 = (a/b)^3, k2 = (a/b)^2. Both reduce exactly to SF = 1 at
      mu_r = 1. The high-permeability asymptotes (WITHOUT the "+1") are
      drawn as dashed lines, reactivating the block that was commented
      out in the original.
  (3) Material annotations corrected: air and polymers have mu_r ~ 1
      (not 10); the annotated points are now mu_r = 1 (non-magnetic,
      SF = 1), 100 (ferrites), 1000 (silicon steel) and 2e4 (Ni-Fe
      alloys / mu-metal).
  (4) Reference fixed: Hoburg (1995) is IEEE Trans. Electromagn. Compat.
      37(4), 574-579, doi 10.1109/15.477342 -- not IEEE Trans. Magn.
  (5) SF is dimensionless and independent of H0: the functions no longer
      depend on the global H0. A horizontal reference line SF = 1 marks
      the no-shielding limit. Unused d removed; results print enabled
      with the correct sphere/cylinder ratio (~1.04 at mu_r = 1000; the
      theoretical maximum of this ratio is 4/3, not 2.2).
"""

import numpy as np
import matplotlib.pyplot as plt

# ===========================================================================================
# Parameters
# ===========================================================================================

mu_r_default = 1000.0  # Default relative permeability for single-point analysis
a = 1.0                # Inner radius (m)
b = 2.0                # Outer radius (m)

# =================================================================
# ANALYTICAL SOLUTIONS FOR SHIELDING FACTOR (EXACT)
# =================================================================

def sphere_shielding_factor(mu_r, a, b):
    """
    Exact shielding factor for a spherical shell, SF = H0/H_int.

    Analytical solution derived from magnetostatic boundary conditions.
    Reference: Jackson, Classical Electrodynamics (3rd ed.), Section 5.12.
    Equivalent closed form:
        SF = (2/9)(1-k3) mu_r + (5+4k3)/9 + 2(1-k3)/(9 mu_r),  k3 = (a/b)^3.
    """
    k3 = a**3 / b**3
    denominator = (mu_r + 2) * (2 * mu_r + 1) - 2 * k3 * (mu_r - 1)**2
    SF = denominator / (9 * mu_r)          # = H0 / H_int (exact)
    return SF


def cylinder_shielding_factor(mu_r, a, b):
    """
    Exact shielding factor for a cylindrical shell (transverse field),
    SF = H0/H_int.

    Analytical solution for an infinite cylindrical shell with the field
    perpendicular to its axis.
    Reference: Hoburg, J.F. (1995), IEEE Trans. Electromagn. Compat.
    37(4), 574-579, doi 10.1109/15.477342.
    Equivalent closed form:
        SF = (1/4)(1-k2) mu_r + (1+k2)/2 + (1-k2)/(4 mu_r),  k2 = (a/b)^2.
    """
    k2 = a**2 / b**2
    denominator = (mu_r + 1)**2 * b**2 - (mu_r - 1)**2 * a**2
    SF = denominator / (4 * mu_r * b**2)   # = H0 / H_int (exact)
    return SF

# =================================================================
# PLOTTING FUNCTION - SHIELDING FACTOR ANALYSIS
# =================================================================

def plot_shielding_factor_analysis(a, b, mu_r_min=1, mu_r_max=50000, num_points=500):
    """
    Plot the exact shielding factor vs relative permeability for both
    geometries, with the high-permeability asymptotes as dashed lines.
    """

    # Generate permeability range (logarithmic scale)
    mu_r_values = np.logspace(np.log10(mu_r_min), np.log10(mu_r_max), num_points)

    # Calculate exact shielding factors
    SF_sphere = sphere_shielding_factor(mu_r_values, a, b)
    SF_cylinder = cylinder_shielding_factor(mu_r_values, a, b)

    # Create figure with single plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))

    # Plot exact shielding factor curves
    ax.loglog(mu_r_values, SF_sphere, 'b-', linewidth=3.5,
              label=f'Spherical shell (a={a:.1f} m, b={b:.1f} m), exact')
    ax.loglog(mu_r_values, SF_cylinder, 'g-', linewidth=3.5,
              label='Cylindrical shell (transverse field), exact')

    # Highlight asymptotic limits (mu_r >> 1), WITHOUT the ad hoc "+1"
    mask = mu_r_values > 100
    SF_sphere_asymptote = (2 / 9) * mu_r_values * (1 - (a**3) / (b**3))
    SF_cylinder_asymptote = (1 / 4) * mu_r_values * (1 - (a**2) / (b**2))
    ax.loglog(mu_r_values[mask], SF_sphere_asymptote[mask],
              'b--', linewidth=1.2, alpha=0.8,
              label=r'Asymptote $\frac{2}{9}\mu_r(1-a^3/b^3)$ (sphere)')
    ax.loglog(mu_r_values[mask], SF_cylinder_asymptote[mask],
              'g--', linewidth=1.2, alpha=0.8,
              label=r'Asymptote $\frac{1}{4}\mu_r(1-a^2/b^2)$ (cylinder)')

    # no-shielding reference
    ax.axhline(1.0, color='gray', linestyle=':', linewidth=1.2, alpha=0.8)

    # Configure plot aesthetics
    ax.set_xlabel('Relative Permeability (μᵣ)', fontsize=20, fontweight='bold')
    ax.set_ylabel('Shielding Factor (SF = H₀/Hᵢₙₜ)', fontsize=20, fontweight='bold')
    ax.set_title('Magnetic Shielding Factor vs Relative Permeability',
                 fontsize=18, fontweight='bold', pad=15)

    ax.grid(True, alpha=0.3, linestyle='--', which='both')
    ax.legend(fontsize=13, framealpha=0.9, loc='upper left')

    # Set axis limits (small padding so the curves do not touch the frame)
    ax.set_xlim(mu_r_min, mu_r_max)
    ax.set_ylim(0.7, 1.5 * max(SF_sphere.max(), SF_cylinder.max()))

    # Reference lines at decades of mu_r
    for mu_r_ref in [1, 10, 100, 1000, 10000]:
        if mu_r_min <= mu_r_ref <= mu_r_max:
            ax.axvline(x=mu_r_ref, color='lightgray', linestyle=':',
                       alpha=0.5, linewidth=0.8)

    # Annotate key permeability values (physically corrected)
    annotation_points = [
        (1, "Non-magnetic media\n(air, polymers): SF = 1", (3, 0.85)),
        (100, "Ferrites\n(initial permeability)", (150, 3.5)),
        (1000, "Silicon steel\n(typical shielding)", (1500, 35)),
        (20000, "High-permeability Ni-Fe\nalloys (mu-metal)", (4000, 2500)),
    ]

    for mu_r_point, text, (tx, ty) in annotation_points:
        if mu_r_min <= mu_r_point <= mu_r_max:
            idx = np.argmin(np.abs(mu_r_values - mu_r_point))
            ax.annotate(text,
                        xy=(mu_r_point, SF_sphere[idx]),
                        xytext=(tx, ty),
                        arrowprops=dict(arrowstyle='->', color='blue', alpha=0.6),
                        fontsize=13,
                        bbox=dict(boxstyle='round,pad=0.3',
                                  facecolor='white', alpha=0.8))

    # Add equation references (EXACT closed forms)
    eq_text = (r'$SF_{\rm sph} = \frac{2}{9}(1-k_3)\mu_r + \frac{5+4k_3}{9}'
               r' + \frac{2(1-k_3)}{9\mu_r}$'
               '\n'
               r'$SF_{\rm cyl} = \frac{1}{4}(1-k_2)\mu_r + \frac{1+k_2}{2}'
               r' + \frac{1-k_2}{4\mu_r}$'
               '\n'
               r'$k_3=(a/b)^3,\ \ k_2=(a/b)^2$   (exact)')

    ax.text(0.98, 0.02, eq_text,
            transform=ax.transAxes,
            verticalalignment='bottom',
            horizontalalignment='right',
            fontsize=14,
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.85))

    plt.tight_layout()
    plt.savefig('Figure_2_modelstyle.png', dpi=150, bbox_inches='tight')
    plt.show()

    # Print detailed analysis
    SF_sph_d = sphere_shielding_factor(mu_r_default, a, b)
    SF_cyl_d = cylinder_shielding_factor(mu_r_default, a, b)
    print("\n" + "=" * 70)
    print("SHIELDING FACTOR ANALYSIS (exact)")
    print("=" * 70)
    print(f"Geometric parameters: a = {a:.3f} m, b = {b:.3f} m")
    print(f"\nShielding factors at mu_r = {mu_r_default:.0f}:")
    print(f"  Spherical shell   : SF = {SF_sph_d:.2f}")
    print(f"  Cylindrical shell : SF = {SF_cyl_d:.2f}")
    print(f"  Ratio (sphere/cylinder) = {SF_sph_d / SF_cyl_d:.3f} "
          f"(theoretical maximum of this ratio: 4/3)")
    print(f"\nChecks: SF(mu_r = 1) sphere = "
          f"{sphere_shielding_factor(1.0, a, b):.6f}, "
          f"cylinder = {cylinder_shielding_factor(1.0, a, b):.6f}  (must be 1)")

    return {
        'mu_r_values': mu_r_values,
        'SF_sphere': SF_sphere,
        'SF_cylinder': SF_cylinder,
        'performance_ratio': SF_sph_d / SF_cyl_d,
    }

# =================================================================
# EXECUTE ANALYSIS
# =================================================================

if __name__ == "__main__":
    results = plot_shielding_factor_analysis(
        a, b,
        mu_r_min=1,
        mu_r_max=50000,
        num_points=500,
    )
"""
Magnetic Shielding Factor Analysis

This code analyzes and plots the shielding factor (S) for spherical and cylindrical shells
as a function of relative permeability. The analysis is based on solutions for
magnetostatic shielding in hollow geometries.
"""

import numpy as np
import matplotlib.pyplot as plt

# ===========================================================================================
# Parameters
# ===========================================================================================

mu_r_default = 1000.0  # Default relative permeability for single-point analysis
H0 = 1.0               # Applied magnetic field intensity (A/m)
a = 1.0                # Inner radius (m)
b = 2.0                # Outer radius (m)  
d = b - a              # Shell thickness (m)

# =================================================================
# ANALYTICAL SOLUTIONS FOR SHIELDING FACTOR
# =================================================================

def sphere_shielding_factor(mu_r, a, b):
    """
    Calculate shielding factor for spherical shell.
    
    Analytical solution derived from magnetostatic boundary conditions.
    Reference: Jackson, Classical Electrodynamics (3rd ed.), Section 5.12
    """
    denominator = (mu_r + 2) * (2 * mu_r + 1) - 2 * (a**3 / b**3) * (mu_r - 1)**2
    H_int = (9 * mu_r * H0) / denominator
    #SF = H0 / H_int
    SF = (2/9)*mu_r*(1-a**3 / b**3)+1
    return SF

def cylinder_shielding_factor(mu_r, a, b):
    """
    Calculate shielding factor for cylindrical shell (transverse field).
    
    Analytical solution for infinite cylindrical shell with field perpendicular to axis.
    Reference: Hoburg, J.F. (1995), IEEE Trans. Magn.
    """
    denominator = (mu_r + 1)**2 * b**2 - (mu_r - 1)**2 * a**2
    H_int = (4 * mu_r * b**2) / denominator * H0
    #SF = H0 / H_int
    SF = (1/4)*mu_r*(1-a**2 / b**2)+1
    return SF

# =================================================================
# PLOTTING FUNCTION - SHIELDING FACTOR ANALYSIS
# =================================================================

def plot_shielding_factor_analysis(a, b, H0=1.0, mu_r_min=1, mu_r_max=50000, num_points=500):
    """
    Plot shielding factor vs relative permeability for both geometries.
    
    Parameters:
    -----------
    a, b : float
        Inner and outer radii (m)
    H0 : float
        External field intensity (A/m)
    mu_r_min, mu_r_max : float
        Range of relative permeability values
    num_points : int
        Number of points in the permeability range
    """
    
    # Generate permeability range (logarithmic scale)
    mu_r_values = np.logspace(np.log10(mu_r_min), np.log10(mu_r_max), num_points)
    
    # Calculate shielding factors
    SF_sphere = np.zeros_like(mu_r_values)
    SF_cylinder = np.zeros_like(mu_r_values)
    
    for i, mu_r in enumerate(mu_r_values):
        SF_sphere[i] = sphere_shielding_factor(mu_r, a, b)
        SF_cylinder[i] = cylinder_shielding_factor(mu_r, a, b)
    
    # Create figure with single plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))
    
    # Plot shielding factor curves
    line_sphere, = ax.loglog(mu_r_values, SF_sphere, 'b-', linewidth=3.5, 
                            label=f'Spherical shell (a={a:.1f}m, b={b:.1f}m)')
    line_cylinder, = ax.loglog(mu_r_values, SF_cylinder, 'g-', linewidth=3.5,
                              label=f'Cylindrical shell (transverse field)')
    
   # # Highlight asymptotic limits
   # if mu_r_max > 100:
        # Calculate asymptotic limits for high permeability
   #    SF_sphere_asymptote = (2/9) * mu_r_values * (1 - (a**3)/(b**3))
  #      SF_cylinder_asymptote = (1/4) * mu_r_values * (1 - (a**2)/(b**2))
        
  #      ax.loglog(mu_r_values[mu_r_values > 100], SF_sphere_asymptote[mu_r_values > 100], 
  #               'b--', linewidth=1.0, alpha=0.7, label='Asymptote: SF ∝ μᵣ (sphere)')
  #      ax.loglog(mu_r_values[mu_r_values > 100], SF_cylinder_asymptote[mu_r_values > 100], 
 #                'g--', linewidth=1.0, alpha=0.7, label='Asymptote: SF ∝ μᵣ (cylinder)')
    
    # Configure plot aesthetics
    ax.set_xlabel('Relative Permeability (μᵣ)', fontsize=20, fontweight='bold')
    ax.set_ylabel('Shielding Factor (SF = H₀/Hᵢₙₜ)', fontsize=20, fontweight='bold')
    ax.set_title('Magnetic Shielding Factor vs Relative Permeability', 
                 fontsize=18, fontweight='bold', pad=15)
    
    ax.grid(True, alpha=0.3, linestyle='--', which='both')
    ax.legend(fontsize=20, framealpha=0.9, loc='upper left')
    
    # Set axis limits
    ax.set_xlim(mu_r_min, mu_r_max)
    y_min = min(np.min(SF_sphere), np.min(SF_cylinder))
    y_max = max(np.max(SF_sphere), np.max(SF_cylinder))
    ax.set_ylim(y_min, y_max)
    
    # Add reference lines and annotations
    mu_r_ref_points = [1, 10, 100, 1000, 10000]
    colors = ['gray', 'lightgray', 'lightgray', 'lightgray', 'lightgray']
    
    for mu_r_ref, color in zip(mu_r_ref_points, colors):
        if mu_r_min <= mu_r_ref <= mu_r_max:
            ax.axvline(x=mu_r_ref, color=color, linestyle=':', alpha=0.5, linewidth=0.8)
    
    # Annotate key permeability values
    annotation_points = [
        (10, "Low-μ materials\n(e.g., air, polymers)"),
        (100, "Ferrites\n(initial permeability)"),
        (1000, "Silicon steel\n(typical shielding)"),
        (10000, "High-permeability\nNi-Fe alloys"),
    ]
    
    for mu_r_point, text in annotation_points:
        if mu_r_min <= mu_r_point <= mu_r_max:
            idx = np.argmin(np.abs(mu_r_values - mu_r_point))
            sf_sphere_val = SF_sphere[idx]
            sf_cyl_val = SF_cylinder[idx]
            
            ax.annotate(text,
                       xy=(mu_r_point, sf_sphere_val),
                       xytext=(mu_r_point*1.5, sf_sphere_val*0.7),
                       arrowprops=dict(arrowstyle='->', color='blue', alpha=0.6),
                       fontsize=16,
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    
    # Calculate and display performance metrics
    idx_default = np.argmin(np.abs(mu_r_values - mu_r_default))
    SF_sphere_default = SF_sphere[idx_default]
    SF_cylinder_default = SF_cylinder[idx_default]
    performance_ratio = SF_sphere_default / SF_cylinder_default
    
    # Add performance comparison box
  #  performance_text = (f'Performance Comparison (μᵣ = {mu_r_default:.0f}):\n'
  #                     f'• Spherical SF = {SF_sphere_default:.1f}\n'
  #                     f'• Cylindrical SF = {SF_cylinder_default:.1f}\n'
  #                     f'• Ratio (sphere/cyl) = {performance_ratio:.2f}\n\n'
  #                     f'Geometric Parameters:\n'
  #                     f'a = {a:.1f} m, b = {b:.1f} m\n'
  #                     f'd = {d:.1f} m, d/(2b) = {d/(2*b):.3f}')
  #  
  #  ax.text(0.02, 0.98, performance_text,
  #          transform=ax.transAxes,
  #          verticalalignment='top',
  #          fontsize=9,
  #          bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    # Add equation references
    eq_text = (r'$SF_{sphere} \approx \frac{2}{9}\mu_r\left(1-\frac{a^3}{b^3}\right)+1$ (μᵣ ≫ 1)'
               '\n'
               r'$SF_{cyl} \approx \frac{1}{4}\mu_r\left(1-\frac{a^2}{b^2}\right)+1$ (μᵣ ≫ 1)')
    
    ax.text(0.98, 0.02, eq_text,
            transform=ax.transAxes,
            verticalalignment='bottom',
            horizontalalignment='right',
            fontsize=16,
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    plt.tight_layout()
    plt.show()
    
    # Print detailed analysis
    #print("\n" + "="*70)
    #print("SHIELDING FACTOR ANALYSIS")
    #print("="*70)
    #print(f"Geometric parameters:")
    #print(f"  Inner radius a = {a:.3f} m")
    #print(f"  Outer radius b = {b:.3f} m")
    #print(f"  Shell thickness d = {d:.3f} m")
    #print(f"  Thickness ratio d/(2b) = {d/(2*b):.4f}")
    
   # print(f"\nShielding factors at μᵣ = {mu_r_default:.0f}:")
    #print(f"  Spherical shell: SF = {SF_sphere_default:.2f}")
   # print(f"  Cylindrical shell: SF = {SF_cylinder_default:.2f}")
   # print(f"  Performance advantage (sphere/cylinder): {performance_ratio:.2f}x")
    
    # Calculate asymptotic values for high permeability
    SF_sphere_asym = (2/9) * mu_r_default * (1 - (a**3)/(b**3))
    SF_cylinder_asym = (1/4) * mu_r_default * (1 - (a**2)/(b**2))
    
    #print(f"\nAsymptotic approximations (μᵣ ≫ 1):")
    #print(f"  Spherical shell: SF ≈ {SF_sphere_asym:.2f}")
    #print(f"  Cylindrical shell: SF ≈ {SF_cylinder_asym:.2f}")
    
   # print(f"\nKey observations:")
   # print(f"  1. Shielding factor increases linearly with μᵣ for μᵣ > 10")
   # print(f"  2. Spherical geometry provides superior shielding effectiveness")
   # print(f"  3. Geometric factor [1-(a/b)³] for sphere is typically larger than [1-(a/b)²] for cylinder")
   # print(f"  4. At μᵣ = {mu_r_default:.0f}, field attenuation inside cavity:")
   # print(f"     • Sphere: H_int/H₀ = {1/SF_sphere_default:.2e} ({100/SF_sphere_default:.3f}% of external field)")
   # print(f"     • Cylinder: H_int/H₀ = {1/SF_cylinder_default:.2e} ({100/SF_cylinder_default:.3f}% of external field)")
    
    return {
        'mu_r_values': mu_r_values,
        'SF_sphere': SF_sphere,
        'SF_cylinder': SF_cylinder,
        'performance_ratio': performance_ratio
    }

# =================================================================
# EXECUTE ANALYSIS
# =================================================================

if __name__ == "__main__":
    results = plot_shielding_factor_analysis(
        a, b, 
        H0=H0, 
        mu_r_min=1, 
        mu_r_max=50000, 
        num_points=500
    )
