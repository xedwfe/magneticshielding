"""
shielding.py
============
Exact static field of a permeable spherical shell and of an infinitely
long cylindrical shell in a uniform applied field, and the quantities
derived from it in the article "Magnetic shielding as a context for
teaching magnetostatics in matter: spherical and cylindrical shells with
open code".

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
  shell; the wall of the shell is the region a < r < b.
* d = 3: spherical shell, applied field along z; fields are evaluated in a
  meridional plane, with in-plane coordinates (u, v) = (z, x).
  d = 2: infinitely long cylindrical shell with axis z, applied field
  along x; fields are evaluated in the transverse plane, (u, v) = (x, y).
* r is the distance from the centre of the sphere or from the axis of the
  cylinder (rho in the article), and theta is measured from the direction
  of the applied field; the poles of a surface are its points with
  theta = 0 and pi, its equator the points with theta = pi/2.
* Field strengths are returned in units of H0 when H0 = 1 (the default).

Offline: depends only on numpy.
"""
import numpy as np

# ---- BEGIN CORE (copied verbatim into the notebook) ----
def first_harmonic_coefficients(mu_r, a, b, d, H0=1.0):
    """Coefficients of the exact solution (Supplementary Material I, section 4).

    Potentials, with r the radial coordinate and t the angle theta measured
    from the applied field:
        cavity   phi_in    = -H_int * r * cos(t)
        shell    phi_shell = (-P * r + Q * r**(1 - d)) * cos(t)
        outside  phi_out   = (-H0 * r + R * r**(1 - d)) * cos(t)
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
    centre = r < 1e-12
    r = np.where(centre, 1e-12, r)
    # at the centre itself theta is undefined; any direction gives the same
    # uniform cavity field, so take theta = 0 there
    c = np.where(centre, 1.0, U / r)
    s = np.where(centre, 0.0, V / r)
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


def flux_function(U, V, mu_r, a, b, d, H0=1.0):
    """Flux function of B (units with mu0 = 1), whose level curves are the
    field lines of B in the plane of fields().

    d = 2: the vector potential A_z, with B_u = dA_z/dv and B_v = -dA_z/du;
           the flux of B per unit length between two field lines is the
           difference of their values.
    d = 3: the Stokes stream function psi = |v| A_phi, multiplied by the
           sign of v, so that the two halves of the meridional plane have
           distinct values; the flux of B through a circle of radius |v|
           about the axis is 2*pi*psi.
    Both have the form G(r) * sin(t) * |sin(t)|**(d - 2), with
        cavity  G = H_int r^(d-1)/(d-1)
        shell   G = mu_r [P r^(d-1)/(d-1) + Q/r]
        outside G = H0 r^(d-1)/(d-1) + R/r,
    which is continuous at r = a and r = b because the normal component of
    B is continuous there.
    """
    U = np.asarray(U, dtype=float)
    V = np.asarray(V, dtype=float)
    r = np.hypot(U, V)
    rs = np.where(r < 1e-12, 1.0, r)
    s = np.where(r < 1e-12, 0.0, V / rs)
    h_int, P, Q, R = first_harmonic_coefficients(mu_r, a, b, d, H0)
    p = r ** (d - 1) / (d - 1)
    G = np.where(r < a, h_int * p,
                 np.where(r <= b, mu_r * (P * p + Q / rs), H0 * p + R / rs))
    return G * s * np.abs(s) ** (d - 2)


def field_line(v0, mu_r, a, b, d, r_max, n=300, H0=1.0):
    """Field line of B through the point (u, v) = (0, v0), v0 != 0, of the
    axis perpendicular to the applied field, in closed form.

    Every field line crosses that axis exactly once, at right angles, so the
    line is the whole level curve of flux_function through (0, v0): with
    F0 = flux_function(0, v0) = G(|v0|) sign(v0), each radius r >= |v0|
    gives |sin(t)|**(d - 1) = |F0|/G(r), that is, two points of the line
    (u = -r cos t upstream and u = +r cos t downstream). No grid and no
    numerical integration are involved, so the line has no gaps and does
    not depend on the resolution of a map. Returns (u, v), ordered along
    the direction of B, for |v0| <= r <= r_max.
    """
    r0 = abs(v0)
    G0 = abs(float(flux_function(0.0, r0, mu_r, a, b, d, H0)))

    def branch(r):
        G = flux_function(np.zeros_like(r), r, mu_r, a, b, d, H0)
        s = np.clip(G0 / G, 0.0, 1.0) ** (1.0 / (d - 1))
        return r * np.sqrt(1.0 - s ** 2), np.sign(v0) * r * s

    w = np.linspace(0.0, 1.0, n)
    r = np.unique(np.concatenate([r0 + (r_max - r0) * w ** 2,
                                  [x for x in (a, b) if r0 < x < r_max]]))
    # refine where consecutive points are far apart: just inside the wall the
    # lines that cross the cavity turn sharply along the inner surface
    for _ in range(40):
        u, v = branch(r)
        big = np.hypot(np.diff(u), np.diff(v)) > r_max / 400
        if not big.any():
            break
        r = np.sort(np.concatenate([r, 0.5 * (r[:-1][big] + r[1:][big])]))
    return np.concatenate([-u[::-1], u[1:]]), np.concatenate([v[::-1], v[1:]])


def field_line_seeds(L, spacing):
    """Points v0 = +-spacing/2, +-3*spacing/2, ... (|v0| < L) of the axis
    perpendicular to the applied field through which the maps draw their
    field lines. Their spacing is a choice of the plot and carries no
    information about the magnitude of the field."""
    n = int(np.floor(L / spacing + 0.5))
    return spacing * (np.arange(-n, n) + 0.5)


def fields_sphere(Z, X, mu_r, a, b, H0=1.0):
    """Spherical shell, meridional plane: returns (Hz, Hx, mu_map, SF)."""
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
    """SF = A*mu_r + B + C/mu_r (article equation (4)); returns (A, B, C),
    with A = C = C_d*(1 - k_d), B = 1 - 2A, C_d = (d-1)/d^2, k_d = x^d."""
    A = (d - 1) / d ** 2 * (1 - x ** d)
    return A, 1 - 2 * A, A


def plus_one_approximation(mu_r, x, d):
    """The common high-permeability approximation SF ~ A*mu_r + 1. It
    exceeds the exact value by A*(2 - 1/mu_r)."""
    A = (d - 1) / d ** 2 * (1 - x ** d)
    return A * mu_r + 1


def regime_crossover(x, d):
    """Permeability near which the plateau SF ~ 1 gives way to SF ~ mu_r,
    mu_r* = 1/[C_d (1 - k_d)] (36/7 for the sphere and 16/3 for the
    cylinder when x = 1/2)."""
    return d ** 2 / ((d - 1) * (1 - x ** d))


def ratio_envelope(x):
    """(SF_sph - 1)/(SF_cyl - 1) = (8/9)(1 + x + x^2)/(1 + x), any mu_r."""
    return 8 / 9 * (1 + x + x ** 2) / (1 + x)


def demagnetizing_factor(d):
    """Demagnetizing factor N = 1/d of the corresponding solid body:
    1/3 for a sphere (d = 3), 1/2 for a long cylinder in a transverse
    field (d = 2). C_d = (d-1)/d^2 = N(1 - N)."""
    return 1.0 / d


demagnetising_factor = demagnetizing_factor   # former name, kept for compatibility


def harmonic_shielding_factor(mu_r, x, m, d):
    """Shielding factor S_m of an applied harmonic of order m
    (Supplementary Material I, section 5): the ratio of its coefficient in
    the applied potential to its coefficient in the cavity, with
    s = m + d - 2. m = 1 is the uniform field; the formula holds for
    m >= 1 (cylinder) and m >= 0 (sphere)."""
    s = m + d - 2
    return 1 + m * s / (m + s) ** 2 * (1 - x ** (m + s)) * (mu_r - 1) ** 2 / mu_r


def cloak_permeability(a, b, d):
    """Relative permeability of a shell of outer radius b, lined by a
    superconducting layer of outer radius a, for which the field outside
    the shell equals the applied field (magnetic cloak, Supplementary
    Material I, section 7, exercise 2)."""
    return ((d - 1) * b ** d + a ** d) / ((d - 1) * (b ** d - a ** d))
