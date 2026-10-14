"""
Is the residual field inside the shielded cavity curved or uniform?
===================================================================

Claim: within the ideal model (linear, homogeneous, isotropic shell;
uniform applied field; full sphere / infinite cylinder), the field in the
cavity is EXACTLY uniform, H = H_int x_hat. The field lines inside are
therefore straight, horizontal, equally spaced. A *small* residual field
does not imply a *curved* one: magnitude and direction are independent
properties.

Four complementary numerical demonstrations:

(A) Harmonic-uniqueness argument.
    The interior solution, regular at the origin, is a_l r^l P_l(cos th)
    (sphere) or a_nu s^nu cos(nu th) (cylinder; s is the cylindrical
    radius, rho in the article, and th is measured from the applied field).
    A UNIFORM applied field contains only the l = 1 (nu = 1) harmonic. For
    every other harmonic, l = 0, 2, 3, ... for the sphere and nu = 0, 2,
    3, ... for the cylinder (nu = 0 with the solutions 1 and ln s), the
    four boundary conditions form a HOMOGENEOUS 4x4 linear system; a
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
    f(R_far) = -H0 R_far on a distant circle). H_int converges at second
    order to the exact solution of the same truncated problem (a 5x5
    linear system). Inside the cavity the computed f(s)/s is constant to
    machine precision, as it must be, since the scheme is exact for f ~ s
    in a homogeneous region.
    This test validates H_int; it cannot test the uniformity, because the
    angular dependence cos(th) is imposed from the start.

(C) Independent two-dimensional solve on a Cartesian grid (transverse
    cylinder), which assumes nothing about the angular dependence.
    div(mu grad phi) = 0 is discretized with the five-point conservative
    scheme, mu taken at the midpoint of each cell face, on a square grid
    of spacing h; the nodes with s >= R_far carry the applied potential
    -H0 x, as in (B). The circular interfaces become staircases, which
    excite spurious higher harmonics of order h, so the computed cavity
    field is not exactly uniform at finite h; its deviation from
    uniformity must fall as the grid is refined, and H_int must approach
    the exact value of the same truncated problem. Requires scipy.

(D) The same test for the spherical shell, on an axisymmetric grid in
    the cylindrical coordinates (s, z) of the meridional half-plane:
    (1/s) d/ds(s mu dphi/ds) + d/dz(mu dphi/dz) = 0 by finite volumes,
    with zero flux through the axis and phi = -H0 z for r >= R_far. No
    angular dependence in theta is assumed; the spherical interfaces
    become staircases. Requires scipy.

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


def det_cylinder0(mu_r, a, b):
    """nu = 0: cavity C, shell E + F ln s, outside D ln s."""
    M = np.array([
        [1.0, -1.0, -np.log(a),  0.0],          # phi continuous at a
        [0.0,  0.0, -mu_r / a,   0.0],          # mu dphi/ds continuous at a
        [0.0,  1.0,  np.log(b), -np.log(b)],    # phi continuous at b
        [0.0,  0.0,  mu_r / b,  -1.0 / b],      # mu dphi/ds continuous at b
    ])
    M = M / np.abs(M).max(axis=1, keepdims=True)
    return np.linalg.det(M)


print(f"{'harmonic':>9s} {'|det| cylinder':>16s} {'|det| sphere':>16s}")
for k in (0, 2, 3, 4, 5, 6):
    dc = det_cylinder0(mu_r, a, b) if k == 0 else det_cylinder(k, mu_r, a, b)
    print(f"{k:>9d} {abs(dc):>16.4e} {abs(det_sphere(k, mu_r, a, b)):>16.4e}")
print("-> all determinants nonzero: only the uniform (l = nu = 1) harmonic")
print("   is excited; the cavity field is EXACTLY uniform.\n")

# =====================================================================
# exact solution of the truncated problem (uniform Dirichlet at R_far)
# unknowns u = (H_i, B2, C2, E3, D3):
#   phi_in = -H_i s cos th ; phi_shell = (B2 s + C2/s) cos th ;
#   phi_out = (E3 s + D3/s) cos th ; phi_out(R_far) = -H0 R_far cos th
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
print("-> H_int converges at second order to the exact value of the truncated")
print("   problem. f(s)/s is constant inside the cavity to machine precision, but")
print("   this is a consistency check, not a test of uniformity: the angular")
print("   dependence cos(th) is imposed here, and the scheme is exact for f ~ s")
print("   in a homogeneous region. Part (C) tests the uniformity without that")
print("   assumption; part (A) shows why the continuous field is exactly uniform.")


# =====================================================================
# (C) two-dimensional Cartesian solve, no assumed angular dependence
# =====================================================================
print()
print("=" * 72)
print("(C) Two-dimensional solve on a Cartesian grid (no assumed cos(th))")
print("=" * 72)
try:
    import scipy.sparse as sps
    import scipy.sparse.linalg as spla
except ImportError:                                  # optional dependency
    sps = None
    print("scipy is not installed: part (C) skipped (pip install scipy).")


def solve_cartesian(h):
    """Five-point conservative scheme on the square [-R_far, R_far]^2;
    unknowns at the nodes with s < R_far, Dirichlet -H0 x elsewhere.
    Returns the mean cavity field H_x and the largest deviation of the
    cavity field vector from that mean, relative to it, both evaluated by
    centred differences at the nodes with s < 0.8 a."""
    n = int(round(2 * R_far / h)) + 1
    g = np.linspace(-R_far, R_far, n)
    X, Y = np.meshgrid(g, g, indexing="ij")
    S = np.hypot(X, Y)
    free = S < R_far - 1e-9
    free[0, :] = free[-1, :] = free[:, 0] = free[:, -1] = False
    phi_bc = -H0 * X

    def mu_at(x, y):
        s_ = np.hypot(x, y)
        return np.where((s_ > a) & (s_ < b), mu_r, 1.0)

    mx = mu_at(0.5 * (X[1:, :] + X[:-1, :]), Y[1:, :])        # faces i+1/2
    my = mu_at(X[:, 1:], 0.5 * (Y[:, 1:] + Y[:, :-1]))        # faces j+1/2
    idx = -np.ones((n, n), dtype=np.int64)
    I, J = np.nonzero(free)
    idx[I, J] = np.arange(I.size)
    cE, cW, cN, cS = mx[I, J], mx[I - 1, J], my[I, J], my[I, J - 1]
    rows, cols, vals = [idx[I, J]], [idx[I, J]], [-(cE + cW + cN + cS)]
    rhs = np.zeros(I.size)
    for di, dj, cc in ((1, 0, cE), (-1, 0, cW), (0, 1, cN), (0, -1, cS)):
        k = idx[I + di, J + dj]
        inside = k >= 0
        rows.append(idx[I, J][inside]); cols.append(k[inside]); vals.append(cc[inside])
        np.add.at(rhs, idx[I, J][~inside], -cc[~inside] * phi_bc[I + di, J + dj][~inside])
    M = sps.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                       shape=(I.size, I.size))
    phi = phi_bc.copy()
    phi[I, J] = spla.spsolve(M.tocsc(), rhs)
    Hx = -(phi[2:, 1:-1] - phi[:-2, 1:-1]) / (2 * h)
    Hy = -(phi[1:-1, 2:] - phi[1:-1, :-2]) / (2 * h)
    core = S[1:-1, 1:-1] < 0.8 * a
    hx, hy = Hx[core], Hy[core]
    mean = hx.mean()
    return mean, np.max(np.hypot(hx - mean, hy)) / abs(mean)


if sps is not None:
    print(f"{'h (m)':>9s} {'h/b':>7s} {'H_int (2D)':>14s} {'vs truncated exact':>19s} "
          f"{'max non-uniformity':>19s}")
    resC = []
    for h in (b / 50, b / 100, b / 200):
        H2, dev = solve_cartesian(h)
        resC.append((h, H2, dev))
        print(f"{h:>9.4f} {h / b:>7.4f} {H2:>14.8e} "
              f"{abs(H2 - H_int_trunc) / H_int_trunc:>19.2e} {dev:>19.2e}")
    print("\n-> with no assumption about the angular dependence, the deviation of the")
    print("   cavity field from uniformity falls steadily under grid refinement")
    print(f"   ({100 * resC[0][2]:.1f}%, {100 * resC[1][2]:.1f}%, {100 * resC[2][2]:.2f}% for h = b/50, b/100, b/200),")
    print("   as H_int approaches the exact value of the truncated problem: the")
    print("   non-uniformity is a discretization error of the staircased interfaces,")
    print("   and the field of the continuous problem is uniform, as (A) shows.")


# =====================================================================
# (D) spherical shell: axisymmetric grid (s, z), no assumed angular law
# =====================================================================
# exact cavity field of the spherical problem truncated at R_far, where
# phi = -H0 z; unknowns (H_i, B2, C2, E3, D3):
#   phi_in = -H_i r cos ; phi_shell = (B2 r + C2/r^2) cos ; phi_out = (E3 r + D3/r^2) cos
M5s = np.array([
    [-a, -a, -1.0 / a**2, 0.0, 0.0],                       # phi continuous at a
    [-1.0, -mu_r, 2 * mu_r / a**3, 0.0, 0.0],              # mu dphi/dr continuous at a
    [0.0, b, 1.0 / b**2, -b, -1.0 / b**2],                 # phi continuous at b
    [0.0, mu_r, -2 * mu_r / b**3, -1.0, 2.0 / b**3],       # mu dphi/dr continuous at b
    [0.0, 0.0, 0.0, R_far, 1.0 / R_far**2],                # Dirichlet at R_far
])
H_int_trunc_sph = np.linalg.solve(M5s, np.array([0.0, 0.0, 0.0, 0.0, -H0 * R_far]))[0]


def solve_axisymmetric(h):
    """Finite volumes on cells of side h in the half-plane s > 0 of the
    meridional plane; cells with r >= R_far carry -H0 z. Returns the mean
    cavity field H_z and the largest relative deviation of the cavity
    field vector from it, at the cells with r < 0.8 a."""
    Ns = int(round(R_far / h))
    s_c = (np.arange(Ns) + 0.5) * h
    z_c = -R_far + (np.arange(2 * Ns) + 0.5) * h
    S, Z = np.meshgrid(s_c, z_c, indexing="ij")
    R = np.hypot(S, Z)
    free = R < R_far
    phi_bc = -H0 * Z

    def mu_at(ss, zz):
        rr = np.hypot(ss, zz)
        return np.where((rr > a) & (rr < b), mu_r, 1.0)

    idx = -np.ones(S.shape, dtype=np.int64)
    I, J = np.nonzero(free)
    k0 = np.arange(I.size)
    idx[I, J] = k0
    rows, cols, vals = [], [], []
    diag = np.zeros(I.size)
    rhs = np.zeros(I.size)
    for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        In, Jn = I + di, J + dj
        valid = (In >= 0) & (In < S.shape[0]) & (Jn >= 0) & (Jn < S.shape[1])
        if di:      # radial face, at s_face, area ~ s_face h (zero on the axis)
            s_face = S[I, J] + 0.5 * di * h
            c = np.where(valid, mu_at(s_face, Z[I, J]) * s_face, 0.0)
        else:       # axial face, area ~ s h
            c = np.where(valid, mu_at(S[I, J], Z[I, J] + 0.5 * dj * h) * S[I, J], 0.0)
        diag -= c
        In_c, Jn_c = np.clip(In, 0, S.shape[0] - 1), np.clip(Jn, 0, S.shape[1] - 1)
        kn = idx[In_c, Jn_c]
        inside = valid & (kn >= 0)
        rows.append(k0[inside]); cols.append(kn[inside]); vals.append(c[inside])
        out = valid & (kn < 0)
        np.add.at(rhs, k0[out], -c[out] * phi_bc[In_c, Jn_c][out])
    rows.append(k0); cols.append(k0); vals.append(diag)
    M = sps.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                       shape=(I.size, I.size))
    phi = phi_bc.copy()
    phi[I, J] = spla.spsolve(M.tocsc(), rhs)
    Hz = -(phi[1:-1, 2:] - phi[1:-1, :-2]) / (2 * h)
    Hs = -(phi[2:, 1:-1] - phi[:-2, 1:-1]) / (2 * h)
    core = R[1:-1, 1:-1] < 0.8 * a
    hz, hs = Hz[core], Hs[core]
    mean = hz.mean()
    return mean, np.max(np.hypot(hz - mean, hs)) / abs(mean)


print()
print("=" * 72)
print("(D) Spherical shell on an axisymmetric grid (no assumed angular law)")
print("=" * 72)
if sps is not None:
    print(f"exact H_int of the sphere, truncated at R = {R_far} m: {H_int_trunc_sph:.8e}")
    print(f"{'h (m)':>9s} {'h/b':>7s} {'H_int (axisym.)':>16s} {'vs truncated exact':>19s} "
          f"{'max non-uniformity':>19s}")
    resD = []
    for h in (b / 50, b / 100, b / 200):
        Hd, dev = solve_axisymmetric(h)
        resD.append((h, Hd, dev))
        print(f"{h:>9.4f} {h / b:>7.4f} {Hd:>16.8e} "
              f"{abs(Hd - H_int_trunc_sph) / H_int_trunc_sph:>19.2e} {dev:>19.2e}")
    print("\n-> for the sphere too, with no assumption about the angular dependence,")
    print("   the deviation of the cavity field from uniformity falls steadily")
    print(f"   ({100 * resD[0][2]:.1f}%, {100 * resD[1][2]:.2f}%, {100 * resD[2][2]:.2f}% for h = b/50, b/100, b/200),")
    print("   and H_int approaches the exact value of the truncated problem.")
else:
    print("scipy is not installed: part (D) skipped (pip install scipy).")
