"""
Is the residual field inside the shielded cavity curved or uniform?
===================================================================

Claim: within the ideal model (linear, homogeneous, isotropic shell;
uniform applied field; full sphere / infinite cylinder), the field in the
cavity is EXACTLY uniform, H = H_int x_hat. The streamlines inside are
therefore straight, horizontal, equally spaced. A *small* residual field
does not imply a *curved* one: magnitude and direction are independent
properties.

Two complementary numerical demonstrations:

(A) Harmonic-uniqueness argument.
    The interior solution, regular at the origin, is a_l r^l P_l(cos th)
    (sphere) or a_nu s^nu cos(nu th) (cylinder). A UNIFORM applied field
    contains only the l = 1 (nu = 1) harmonic. For every other harmonic
    the four boundary conditions form a HOMOGENEOUS 4x4 linear system; a
    nonzero determinant means that harmonic is not excited at all. The
    determinants below never vanish. Hence the cavity potential can only
    be  -H_int r cos(th):  a uniform field. Curvature of the interior
    lines would require some other harmonic -- and there is none.

(B) Independent radial solve of the nu = 1 mode (transverse cylinder).
    Given (A), the full potential is phi = f(s) cos(th) and the PDE
    div(mu grad phi) = 0 reduces to  (s mu f')' - mu f / s = 0.  We solve
    this two-point boundary-value problem by a conservative finite-volume
    scheme that knows nothing about the analytic solution (interfaces at
    s = a, b fall exactly on cell faces; regularity f(0) = 0; Dirichlet
    f(R_far) = -H0 R_far on a distant circle). Inside the cavity the
    computed f(s)/s is constant to machine precision (f is a straight line
    through the origin => uniform field), the deviation falls as h^2
    under grid refinement, and H_int matches the exact solution of the
    same truncated problem (a 5x5 linear system) to the same accuracy.
    (A 2D Cartesian grid is deliberately avoided: it staircases the
    circular interfaces and injects spurious high harmonics -- a
    discretization artifact, not cavity physics.)

Where WOULD curvature come from in a real shield? Only from ingredients
outside this model: gradients of the applied field (they source l > 1
directly), finite cylinder length / open ends, holes and seams, magnetic
saturation (mu depending on H), inhomogeneity or anisotropy of the
material. Good discussion points for the paper.
"""

import numpy as np

mu_r = 1000.0
a, b = 1.0, 2.0
H0 = 1.0
R_far = 5.5          # radius of the circular far boundary (truncated problem)

# =====================================================================
# (A) determinants of the homogeneous boundary-condition systems
# =====================================================================
print("=" * 72)
print("(A) Harmonic-uniqueness check: |det| of the homogeneous BC system")
print("    (nonzero determinant  =>  that harmonic is NOT excited)")
print("=" * 72)


def det_cylinder(nu, mu_r, a, b):
    M = np.array([
        [a**nu,            -a**nu,             -a**(-nu),               0.0],
        [a**(nu - 1), -mu_r * a**(nu - 1),  mu_r * a**(-nu - 1),        0.0],
        [0.0,               b**nu,              b**(-nu),        -b**(-nu)],
        [0.0,          mu_r * b**(nu - 1), -mu_r * b**(-nu - 1),  b**(-nu - 1)],
    ])
    M = M / np.abs(M).max(axis=1, keepdims=True)
    return np.linalg.det(M)


def det_sphere(l, mu_r, a, b):
    M = np.array([
        [a**l,             -a**l,                    -a**(-(l + 1)),                 0.0],
        [l * a**(l - 1), -mu_r * l * a**(l - 1),  mu_r * (l + 1) * a**(-(l + 2)),    0.0],
        [0.0,               b**l,                     b**(-(l + 1)),         -b**(-(l + 1))],
        [0.0,          mu_r * l * b**(l - 1), -mu_r * (l + 1) * b**(-(l + 2)),
         (l + 1) * b**(-(l + 2))],
    ])
    M = M / np.abs(M).max(axis=1, keepdims=True)
    return np.linalg.det(M)


print(f"{'harmonic':>9s} {'|det| cylinder':>16s} {'|det| sphere':>16s}")
for k in range(2, 7):
    print(f"{k:>9d} {abs(det_cylinder(k, mu_r, a, b)):>16.4e} "
          f"{abs(det_sphere(k, mu_r, a, b)):>16.4e}")
print("-> all determinants nonzero: only the uniform (l = nu = 1) harmonic")
print("   is excited; the cavity field is EXACTLY uniform.\n")

# =====================================================================
# exact solution of the truncated problem (uniform Dirichlet at R_far)
# unknowns u = (H_i, B2, C2, E3, D3):
#   phi1 = -H_i s cos th ; phi2 = (B2 s + C2/s) cos th ;
#   phi3 = (E3 s + D3/s) cos th ; phi3(R_far) = -H0 R_far cos th
# =====================================================================
M5 = np.array([
    [-a, -a, -1.0 / a, 0.0, 0.0],                    # phi continuous at a
    [1.0, mu_r, -mu_r / a**2, 0.0, 0.0],             # B_n continuous at a
    [0.0, b, 1.0 / b, -b, -1.0 / b],                 # phi continuous at b
    [0.0, mu_r, -mu_r / b**2, -1.0, 1.0 / b**2],     # B_n continuous at b
    [0.0, 0.0, 0.0, R_far, 1.0 / R_far],             # Dirichlet at R_far
])
rhs = np.array([0.0, 0.0, 0.0, 0.0, -H0 * R_far])
H_int_trunc = np.linalg.solve(M5, rhs)[0]

denom = (mu_r + 1) ** 2 * b**2 - (mu_r - 1) ** 2 * a**2
H_int_inf = 4 * mu_r * b**2 * H0 / denom

# =====================================================================
# (B) conservative finite-volume solve of (s mu f')' - mu f / s = 0
# =====================================================================
print("=" * 72)
print("(B) Independent radial solve of the nu = 1 mode (finite volumes)")
print("=" * 72)
print(f"exact H_int, infinite domain        : {H_int_inf:.8e}")
print(f"exact H_int, truncated at R = {R_far} m : {H_int_trunc:.8e}   "
      f"(finite-boundary effect: {abs(H_int_trunc - H_int_inf) / H_int_inf:.1%})")
print()


def solve_radial(h):
    """Cell-centered FV grid on (0, R_far]; a, b and R_far are exact
    multiples of h, so every material interface lies on a cell face."""
    Ncell = int(round(R_far / h))
    faces = h * np.arange(Ncell + 1)            # s = 0 ... R_far
    centers = 0.5 * (faces[:-1] + faces[1:])

    mu_cell = np.ones(Ncell)
    mu_cell[(centers > a) & (centers < b)] = mu_r
    # face permeabilities: harmonic mean of the neighbouring cells
    mu_face = np.empty(Ncell + 1)
    mu_face[1:-1] = 2 * mu_cell[:-1] * mu_cell[1:] / (mu_cell[:-1] + mu_cell[1:])
    mu_face[0], mu_face[-1] = mu_cell[0], mu_cell[-1]

    # tridiagonal system A f = r for f at cell centers
    lo = np.zeros(Ncell)      # sub-diagonal
    di = np.zeros(Ncell)      # diagonal
    up = np.zeros(Ncell)      # super-diagonal
    r = np.zeros(Ncell)

    wE = faces[1:] * mu_face[1:] / h            # east-face conductances
    wW = faces[:-1] * mu_face[:-1] / h          # west face (zero at s = 0)
    sink = mu_cell * h / centers                # from -mu f / s term

    di[:] = -(wE + wW) - sink
    up[:-1] = wE[:-1]
    lo[1:] = wW[1:]
    # outer Dirichlet f(R_far) = -H0 R_far via ghost value at the last face
    di[-1] -= wE[-1]                            # ghost: f_ghost = 2 f_bc - f_N
    r[-1] = -2.0 * wE[-1] * (-H0 * R_far)

    # Thomas algorithm
    for j in range(1, Ncell):
        m = lo[j] / di[j - 1]
        di[j] -= m * up[j - 1]
        r[j] -= m * r[j - 1]
    f = np.empty(Ncell)
    f[-1] = r[-1] / di[-1]
    for j in range(Ncell - 2, -1, -1):
        f[j] = (r[j] - up[j] * f[j + 1]) / di[j]

    core = centers < 0.9 * a
    g = f[core] / centers[core]                 # = -H_int if f is linear
    H_int_fd = -g.mean()
    lin_dev = np.max(np.abs(g / g.mean() - 1.0))
    return H_int_fd, lin_dev


print(f"{'h (m)':>9s} {'H_int (FD)':>14s} {'vs truncated exact':>19s} "
      f"{'max dev. of f(s)/s':>19s}")
res = []
for h in (0.01, 0.005, 0.0025):
    H_fd, lin_dev = solve_radial(h)
    res.append((h, H_fd, lin_dev))
    print(f"{h:>9.4f} {H_fd:>14.8e} "
          f"{abs(H_fd - H_int_trunc) / H_int_trunc:>19.2e} {lin_dev:>19.2e}")

err = [abs(H - H_int_trunc) / H_int_trunc for _, H, _ in res]
print(f"\nrefinement ratios of the H_int error: {err[0] / err[1]:.2f}, "
      f"{err[1] / err[2]:.2f}  (~4 = clean second-order convergence)")
print("note: the deviation of f(s)/s from a constant is already at MACHINE")
print("PRECISION (~1e-13) at every resolution -- the computed cavity field")
print("is uniform to rounding error, independent of h.")
print("-> f(s) is a straight line through the origin inside the cavity:")
print("   the interior field is UNIFORM and the streamlines are STRAIGHT.")
print("   Curvature would require l > 1 harmonics, which (A) shows are")
print("   not excited by a uniform applied field.")
