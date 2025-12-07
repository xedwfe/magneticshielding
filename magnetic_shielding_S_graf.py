"""
Magnetic Shielding Factor Analysis

This code analyzes and plots the shielding factor (SF) for spherical and cylindrical shells
as a function of relative permeability. The analysis is based on analytical solutions for
magnetostatic shielding in hollow geometries.
"""

import numpy as np
import matplotlib.pyplot as plt

# =================================================================
# Simulation Parameters
# =================================================================

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
    SF = H0 / H_int
    return SF

def cylinder_shielding_factor(mu_r, a, b):
    """
    Calculate shielding factor for cylindrical shell (transverse field).
    
    Analytical solution for infinite cylindrical shell with field perpendicular to axis.
    Reference: Hoburg, J.F. (1995), IEEE Trans. Magn.
    """
    denominator = (mu_r + 1)**2 * b**2 - (mu_r - 1)**2 * a**2
    H_int = (4 * mu_r * b**2) / denominator * H0
    SF = H0 / H_int
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
    line_sphere, = ax.loglog(mu_r_values, SF_sphere, 'b-', linewidth=2.5, 
                            label=f'Spherical shell (a={a:.1f}m, b={b:.1f}m)')
    line_cylinder, = ax.loglog(mu_r_values, SF_cylinder, 'g-', linewidth=2.5,
                              label=f'Cylindrical shell (transverse field)')
    
    # Highlight asymptotic limits
    if mu_r_max > 100:
        # Calculate asymptotic limits for high permeability
        SF_sphere_asymptote = (2/9) * mu_r_values * (1 - (a**3)/(b**3))
        SF_cylinder_asymptote = (1/4) * mu_r_values * (1 - (a**2)/(b**2))
        
        ax.loglog(mu_r_values[mu_r_values > 100], SF_sphere_asymptote[mu_r_values > 100], 
                 'b--', linewidth=1.0, alpha=0.7, label='Asymptote: SF ∝ μᵣ (sphere)')
        ax.loglog(mu_r_values[mu_r_values > 100], SF_cylinder_asymptote[mu_r_values > 100], 
                 'g--', linewidth=1.0, alpha=0.7, label='Asymptote: SF ∝ μᵣ (cylinder)')
    
    # Configure plot aesthetics
    ax.set_xlabel('Relative Permeability (μᵣ)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Shielding Factor (SF = H₀/Hᵢₙₜ)', fontsize=12, fontweight='bold')
    ax.set_title('Magnetic Shielding Factor vs Relative Permeability', 
                 fontsize=14, fontweight='bold', pad=15)
    
    ax.grid(True, alpha=0.3, linestyle='--', which='both')
    ax.legend(fontsize=10, framealpha=0.9, loc='upper left')
    
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
                       fontsize=8,
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
    eq_text = (r'$SF_{sphere} \approx \frac{2}{9}\mu_r\left(1-\frac{a^3}{b^3}\right)$ (μᵣ ≫ 1)'
               '\n'
               r'$SF_{cyl} \approx \frac{1}{4}\mu_r\left(1-\frac{a^2}{b^2}\right)$ (μᵣ ≫ 1)')
    
    ax.text(0.98, 0.02, eq_text,
            transform=ax.transAxes,
            verticalalignment='bottom',
            horizontalalignment='right',
            fontsize=9,
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