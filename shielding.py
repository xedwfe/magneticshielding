"""
shielding.py
============
Exact static field of a permeable spherical shell and of an infinitely
long cylindrical shell in a uniform applied field, and the quantities
derived from it in the article "Magnetic shielding as a motivational tool
in teaching classical electromagnetic theory".

This module is the single source of the physics used by paper_figures.py
and explore_shielding.py. The interactive notebook carries an identical
copy of the block between the CORE markers, so that it runs on Google
Colab without any other file; verify_revision_claims.py checks that the
two copies agree and compares this module with an independent
implementation.

Model and conventions (as in the article)
-----------------------------------------
* Linear, isotropic and homogeneous shell of relative permeability mu_r,
  inner radius a and outer radius b, in vacuum, in a static uniform
  applied field H0; far from the shell B0 = mu0*H0.
* B = mu0*H in the cavity and outside the shell, B = mu0*mu_r*H in the
  shell.
* d = 3: spherical shell, applied field along z; fields are evaluated in a
  meridian plane, with in-plane coordinates (u, v) = (z, x).
  d = 2: infinitely long cylindrical shell with axis z, applied field
  along x; fields are evaluated in the transverse plane, (u, v) = (x, y).
* Field strengths are returned in units of H0 when H0 = 1 (the default).

Offline: depends only on numpy.
"""
import numpy as np

# ---- BEGIN CORE (copied verbatim into the notebook) ----
def first_harmonic_coefficients(mu_r, a, b, d, H0=1.0):
    """Coefficients of the exact solution (Supplementary Material I, section 4).

    Potentials, with r the distance to the centre and t the angle to the
    applied field:
        cavity   phi1 = -H_int * r * cos(t)
        shell    phi2 = (-P * r + Q * r**(1 - d)) * cos(t)
        outside  phi3 = (-H0 * r + R * r**(1 - d)) * cos(t)
    Returns (H_int, P, Q, R).
    """
    x = a / b
    delta = (mu_r + d - 1) * ((d - 1) * mu_r + 1) - (d - 1) * (mu_r - 1) ** 2 * x ** d
    h_int = d ** 2 * mu_r * H0 / delta
    P = d * ((d - 1) * mu_r + 1) * H0 / delta
    Q = -d * (mu_r - 1) * a ** d * H0 / delta
    R = (mu_r - 1) * ((d - 1) * mu_r + 1) * (b ** d - a ** d) * H0 / delta
    return h_int, P, Q, R


def fields(U, V, mu_r, a, b, d, H0=1.0):
    """Field strength at the points (U, V) of the plane described above.

    Returns (Hu, Hv, mu_map, SF): the components of H along u (the applied
    field direction) and v, the map of the relative permeability (mu_r in
    the shell, 1 elsewhere) and the shielding factor SF = H0/H_int.
    The flux density is B = mu0 * mu_map * H.
    """
    U = np.asarray(U, dtype=float)
    V = np.asarray(V, dtype=float)
    r = np.hypot(U, V)
    r = np.where(r < 1e-12, 1e-12, r)
    c, s = U / r, V / r
    h_int, P, Q, R = first_harmonic_coefficients(mu_r, a, b, d, H0)
    cavity = r < a
    shell = (r >= a) & (r <= b)
    Hr = np.where(cavity, h_int * c,
                  np.where(shell, (P + (d - 1) * Q * r ** -d) * c,
                           (H0 + (d - 1) * R * r ** -d) * c))
    Ht = np.where(cavity, -h_int * s,
                  np.where(shell, (-P + Q * r ** -d) * s,
                           (-H0 + R * r ** -d) * s))
    Hu = Hr * c - Ht * s
    Hv = Hr * s + Ht * c
    mu_map = np.where(shell, mu_r, 1.0)
    return Hu, Hv, mu_map, H0 / h_int


def fields_sphere(Z, X, mu_r, a, b, H0=1.0):
    """Spherical shell, meridian plane: returns (Hz, Hx, mu_map, SF)."""
    return fields(Z, X, mu_r, a, b, 3, H0)


def fields_cylinder(X, Y, mu_r, a, b, H0=1.0):
    """Cylindrical shell, transverse plane: returns (Hx, Hy, mu_map, SF)."""
    return fields(X, Y, mu_r, a, b, 2, H0)


def shielding_factor(mu_r, x, d):
    """Exact SF = H0/H_int, article equation (3), with x = a/b."""
    return 1 + (d - 1) / d ** 2 * (1 - x ** d) * (mu_r - 1) ** 2 / mu_r


# a/b at which the sphere and the cylinder shield equally, for every mu_r
X_STAR = (1 + np.sqrt(33)) / 16
# ---- END CORE ----


def SF_sph(mu_r, x):
    """Exact shielding factor of the spherical shell."""
    return shielding_factor(mu_r, x, 3)


def SF_cyl(mu_r, x):
    """Exact shielding factor of the transverse cylindrical shell."""
    return shielding_factor(mu_r, x, 2)


def three_term_coefficients(x, d):
    """SF = A*mu_r + B + C/mu_r (article equation (4)); returns (A, B, C)."""
    A = (d - 1) / d ** 2 * (1 - x ** d)
    return A, 1 - 2 * A, A


def plus_one_approximation(mu_r, x, d):
    """The common high-permeability approximation SF ~ A*mu_r + 1. It
    exceeds the exact value by A*(2 - 1/mu_r)."""
    A = (d - 1) / d ** 2 * (1 - x ** d)
    return A * mu_r + 1


def regime_crossover(x, d):
    """Permeability near which the plateau SF ~ 1 gives way to SF ~ mu_r."""
    return d ** 2 / ((d - 1) * (1 - x ** d))


def ratio_envelope(x):
    """(SF_sph - 1)/(SF_cyl - 1) = (8/9)(1 + x + x^2)/(1 + x), any mu_r."""
    return 8 / 9 * (1 + x + x ** 2) / (1 + x)


def demagnetising_factor(d):
    """N = 1/3 for a sphere (d = 3), 1/2 for a long cylinder in a
    transverse field (d = 2)."""
    return 1.0 / d


def harmonic_shielding_factor(mu_r, x, m, d):
    """Shielding factor of an applied harmonic of order m (Supplementary
    Material I, section 5); m = 1 is the uniform field."""
    s = m + d - 2
    return 1 + m * s / (m + s) ** 2 * (1 - x ** (m + s)) * (mu_r - 1) ** 2 / mu_r


def cloak_permeability(R1, R2, d):
    """Permeability of a shell (radii R1 < R2) around a superconducting core
    that leaves the external field undisturbed (Supplementary Material I,
    section 7)."""
    return ((d - 1) * R2 ** d + R1 ** d) / ((d - 1) * (R2 ** d - R1 ** d))
