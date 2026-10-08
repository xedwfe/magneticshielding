"""Offline numerical verification of every quantitative statement made in the
revised manuscript EJP-110956 and its supplements, and of the consistency of
the repository code (shielding.py, paper_figures.py and the notebook).
Run from the repository folder:  python verify_revision_claims.py
No network access is used."""
import numpy as np
ok_all = True
def check(name, cond):
    global ok_all; ok_all &= bool(cond); print(f"[{'OK' if cond else 'FAIL'}] {name}")

def SF_sph(mu, x): k = x**3; return (2/9)*(1-k)*mu + (5+4*k)/9 + 2*(1-k)/(9*mu)
def SF_cyl(mu, x): k = x**2; return (1/4)*(1-k)*mu + (1+k)/2 + (1-k)/(4*mu)
mu = np.logspace(-3, 5, 600); xs = np.linspace(0.02, 0.98, 49)

# 1. unified compact form SF = 1 + (d-1)/d^2 [1-x^d] (mu-1)^2/mu
for d, f in [(3, SF_sph), (2, SF_cyl)]:
    check(f"compact form, d={d}", all(np.allclose(f(mu, x), 1+(d-1)/d**2*(1-x**d)*(mu-1)**2/mu) for x in xs))

# 2. exact solve of the 4x4 systems (cylinder nu, sphere l) vs closed forms, incl. higher harmonics
def S_cyl_nu(mu, a, b, n):
    # unknowns C (int), E, F (shell), D (ext); applied -G rho^n cos(n th), G=1
    M = np.array([[a**n, -a**n, -a**-n, 0],
                  [a**n, -mu*a**n, mu*a**-n, 0],
                  [0, b**n, b**-n, -b**-n],
                  [0, mu*b**n, -mu*b**-n, b**-n]], float)
    rhs = np.array([0, 0, -b**n, -b**n], float)
    C = np.linalg.solve(M, rhs)[0]; return 1/abs(C)
def S_sph_l(mu, a, b, l):
    M = np.array([[a**l, -a**l, -a**-(l+1), 0],
                  [l*a**l, -mu*l*a**l, mu*(l+1)*a**-(l+1), 0],
                  [0, b**l, b**-(l+1), -b**-(l+1)],
                  [0, mu*l*b**l, -mu*(l+1)*b**-(l+1), (l+1)*b**-(l+1)]], float)
    rhs = np.array([0, 0, -b**l, -l*b**l], float)
    C = np.linalg.solve(M, rhs)[0]; return 1/abs(C)
a, b = 1.0, 2.0
for m in [0.01, 1, 3.7, 10, 1e3]:
    check(f"4x4 solve = SF_cyl, mu={m}", np.isclose(S_cyl_nu(m, a, b, 1), SF_cyl(m, a/b)))
    check(f"4x4 solve = SF_sph, mu={m}", np.isclose(S_sph_l(m, a, b, 1), SF_sph(m, a/b)))
    for n in [2, 3, 5]:
        x = a/b
        Sc = ((m+1)**2 - (m-1)**2*x**(2*n))/(4*m)
        Ss = ((l:=n)*m + l + 1)*((l+1)*m + l) - l*(l+1)*(m-1)**2*x**(2*l+1)
        Ss = Ss/((2*l+1)**2*m)
        check(f"  harmonic {n}: closed forms (cyl, sph), mu={m}",
              np.isclose(S_cyl_nu(m, a, b, n), Sc) and np.isclose(S_sph_l(m, a, b, n), Ss))
# higher harmonics shielded better (Bidinosti & Martin 2014)
m = 1000.
check("higher harmonics shielded progressively better",
      all(S_cyl_nu(m,a,b,n+1) > S_cyl_nu(m,a,b,n) and S_sph_l(m,a,b,n+1) > S_sph_l(m,a,b,n) for n in range(1,6)))

# 3. homogeneous systems: determinants never vanish
def det_c(m, x, n): return (m+1)**2 - (m-1)**2*x**(2*n)
def det_s(m, x, l): return (l*m+l+1)*((l+1)*m+l) - l*(l+1)*(m-1)**2*x**(2*l+1)
grid = [(m, x, n) for m in np.logspace(-4, 6, 60) for x in np.linspace(0.001, 0.999, 60) for n in range(1, 9)]
check("cyl determinant > 0 for all mu>0, 0<a<b, nu>=1", min(det_c(*g) for g in grid) > 0)
check("sph determinant > 0 for all mu>0, 0<a<b, l>=1", min(det_s(*g) for g in grid) > 0)

# 4. geometry comparison
xstar = (1+np.sqrt(33))/16
check(f"x* = (1+sqrt33)/16 = {xstar:.4f} solves 8x^2-x-1=0", abs(8*xstar**2-xstar-1) < 1e-14)
for m in [1.5, 10, 100, 1e3, 1e5]:
    check(f"SF_sph = SF_cyl exactly at x*, mu={m}", np.isclose(SF_sph(m, xstar), SF_cyl(m, xstar)))
r = (SF_sph(37., xs)-1)/(SF_cyl(37., xs)-1)
check("(SF_sph-1)/(SF_cyl-1) = (8/9)(1+x+x^2)/(1+x), mu-independent", np.allclose(r, (8/9)*(1+xs+xs**2)/(1+xs)))
R = SF_sph(mu[:,None], xs[None,:])/SF_cyl(mu[:,None], xs[None,:])
env = (8/9)*(1+xs+xs**2)/(1+xs)
check("SF_sph/SF_cyl lies between 1 and the envelope (8/9..4/3)",
      np.all((R-1)*(R-env) <= 1e-12))
print(f"      a/b=0.30, mu=1000: SF_sph={SF_sph(1000,.3):.1f}  SF_cyl={SF_cyl(1000,.3):.1f}")
print(f"      a/b=0.50, mu=1000: SF_sph={SF_sph(1000,.5):.2f} SF_cyl={SF_cyl(1000,.5):.2f} ratio={SF_sph(1000,.5)/SF_cyl(1000,.5):.3f}")
print(f"      a/b=0.80, mu=1000: SF_sph={SF_sph(1000,.8):.1f}  SF_cyl={SF_cyl(1000,.8):.1f}")

# 5. demagnetising-factor reading, N = 1/d
for d, f in [(3, SF_sph), (2, SF_cyl)]:
    N = 1/d; m = 1e7
    check(f"thin-shell: (SF-1)/(mu t/b) -> 1-N, d={d}", np.isclose((f(m, 1-1e-5)-1)/(m*1e-5), 1-N, rtol=1e-3))
    check(f"thick-shell: (SF-1)/mu -> N(1-N), d={d}", np.isclose((f(m, 1e-4)-1)/m, N*(1-N), rtol=1e-3))
    Hsolid = 1/(1+N*(m-1)); Hcav = Hsolid*m/(m*(1-N)+N)
    check(f"solid body x cavity amplification gives N(1-N)mu, d={d}", np.isclose(1/Hcav, N*(1-N)*m, rtol=1e-3))

# 6. '+1' approximation: excess A(2-1/mu); relative values for a/b = 1/2
A = (2/9)*(1-0.125); mm = np.linspace(1, 10, 900001)
rel = (A*mm+1)/SF_sph(mm, .5) - 1
check("absolute excess = A(2-1/mu)", np.allclose(A*mm+1-SF_sph(mm,.5), A*(2-1/mm)))
print(f"      relative excess: mu=1 {100*rel[0]:.1f}%, max {100*rel.max():.1f}% at mu={mm[rel.argmax()]:.2f}, mu=10 {100*rel[-1]:.1f}%")
check("relative excess ~ 2/mu at large mu", np.isclose((A*1e4+1)/SF_sph(1e4,.5)-1, 2e-4, rtol=1e-3))

# 7. quantified attenuation for a/b = 1/2 (replaces the 'order 10^3' threshold)
for m in [10, 100, 1000]:
    print(f"      mu={m:>5}: SF_sph={SF_sph(m,.5):7.2f}  H_int/H0={1/SF_sph(m,.5):.4f}   SF_cyl={SF_cyl(m,.5):7.2f}")
mu_star_s, mu_star_c = 9/(2*(1-.125)), 4/(1-.25)
print(f"      crossover mu*: sphere {mu_star_s:.2f}, cylinder {mu_star_c:.2f}")

# 8. mu <-> 1/mu symmetry and superconducting limit
check("SF(mu) = SF(1/mu) for both geometries", np.allclose(SF_sph(mu,.5), SF_sph(1/mu,.5)) and np.allclose(SF_cyl(mu,.5), SF_cyl(1/mu,.5)))

# 9. magnetic cloak: superconducting core (R1) + permeable shell (R1..R2), no exterior perturbation
def cloak_D(mu, R1, R2, geom):
    if geom == "cyl":   # shell (E rho + F/rho)cos, ext (-rho + D/rho)cos ; B_n=0 at R1
        M = np.array([[1, -1/R1**2, 0], [R2, 1/R2, -1/R2], [mu, -mu/R2**2, 1/R2**2]])
        rhs = np.array([0, -R2, -1.]); return np.linalg.solve(M, rhs)[2]
    else:               # shell (E r + F/r^2)cos, ext (-r + D/r^2)cos
        M = np.array([[1, -2/R1**3, 0], [R2, 1/R2**2, -1/R2**2], [mu, -2*mu/R2**3, 2/R2**3]])
        rhs = np.array([0, -R2, -1.]); return np.linalg.solve(M, rhs)[2]
R1, R2 = 1.0, 1.3
mc = (R2**2+R1**2)/(R2**2-R1**2); ms = (2*R2**3+R1**3)/(2*(R2**3-R1**3))
check(f"cylinder cloak: mu=(R2^2+R1^2)/(R2^2-R1^2)={mc:.3f} -> no exterior dipole", abs(cloak_D(mc,R1,R2,'cyl')) < 1e-12)
check(f"sphere cloak:  mu=(2R2^3+R1^3)/(2(R2^3-R1^3))={ms:.3f} -> no exterior dipole", abs(cloak_D(ms,R1,R2,'sph')) < 1e-12)


# 10. single derivation for both geometries with d as a parameter (Supplementary Material I)
#     interior  -H_int r cos ; shell (-P r + Q r^(1-d)) cos ; exterior (-H0 r + R r^(1-d)) cos
ok10 = True
for d in (2, 3):
    for m in (0.3, 1.0, 7.0, 1e3):
        for x in (0.1, 0.5, 0.9):
            a, b, H0 = x, 1.0, 1.0
            M = np.array([[1, -1, a**-d, 0],                 # phi continuous at a:  H_int = P - Q a^-d
                          [1, -m, -m*(d-1)*a**-d, 0],         # mu dphi/dr at a:     H_int = mu(P + (d-1)Q a^-d)
                          [0, 1, -b**-d, b**-d],              # phi continuous at b:  P - Q b^-d = H0 - R b^-d
                          [0, m, m*(d-1)*b**-d, -(d-1)*b**-d]]) # mu dphi/dr at b:   mu(P+(d-1)Q b^-d) = H0+(d-1)R b^-d
            Hint, P, Q, R = np.linalg.solve(M, np.array([0, 0, H0, H0]))
            Dd = (m + d - 1)*((d - 1)*m + 1) - (d - 1)*(m - 1)**2*x**d
            ok10 &= np.isclose(Hint, d*d*m*H0/Dd) and np.isclose(P, d*((d-1)*m+1)*H0/Dd) \
                and np.isclose(Q, d*(1-m)*a**d*H0/Dd) and np.isclose(R, (m-1)*((d-1)*m+1)*(b**d-a**d)*H0/Dd)
            ok10 &= np.isclose(Dd, d*d*m + (d-1)*(m-1)**2*(1-x**d))
check("general-d coefficients H_int, P, Q, R and identity for Delta_d", ok10)

# 11. unified formula for any driven harmonic m (cylinder: s = m; sphere: s = m+1)
ok11 = True
for m_ in (0.05, 1.0, 4.0, 300.0):
    for x in (0.2, 0.5, 0.9):
        for n in range(1, 7):
            for d, S in ((2, S_cyl_nu), (3, S_sph_l)):
                mm, s = n, n + d - 2
                uni = 1 + mm*s/(mm+s)**2*(1-x**(mm+s))*(m_-1)**2/m_
                ok11 &= np.isclose(S(m_, x, 1.0, n), uni)
                det = (mm+m_*s)*(s+m_*mm) - mm*s*(m_-1)**2*x**(mm+s)
                ok11 &= det >= m_*(mm+s)**2*(1-1e-12)
check("S_m = 1 + ms/(m+s)^2 (1-x^(m+s))(mu-1)^2/mu and Det >= mu (m+s)^2 b^(m+s)", ok11)

# 12. exact small-cavity limit from demagnetising factors, N = 1/d
for d in (2, 3):
    N = 1/d; mus = np.logspace(-2, 4, 50)
    check(f"SF(a->0) = [1+N(mu-1)][1+(1-N)(mu-1)]/mu exactly, d={d}",
          np.allclose(1 + (d-1)/d**2*(mus-1)**2/mus, (1+N*(mus-1))*(1+(1-N)*(mus-1))/mus))

# 13. unified cloak condition mu = [(d-1)R2^d + R1^d]/[(d-1)(R2^d - R1^d)]
for d, g in ((2, "cyl"), (3, "sph")):
    R1, R2 = 1.0, 1.3
    mu_cl = ((d-1)*R2**d + R1**d)/((d-1)*(R2**d - R1**d))
    check(f"unified cloak condition, d={d}: mu={mu_cl:.3f}", abs(cloak_D(mu_cl, R1, R2, g)) < 1e-12)

# 14. the fields actually plotted (paper_figures.py) satisfy the interface conditions
import paper_figures as pf
okI = True
for geom, f in (("sph", pf.fields_sphere), ("cyl", pf.fields_cylinder)):
    for mu_ in (3.0, 1000.0):
        for R0 in (1.0, 2.0):
            th = np.linspace(0.1, 3.0, 13)
            def comp(R):
                U, V = R*np.cos(th), R*np.sin(th)
                Hu, Hv, mr, _ = f(U, V, mu_, 1.0, 2.0)
                n = np.array([np.cos(th), np.sin(th)]); t = np.array([-np.sin(th), np.cos(th)])
                return Hu*t[0] + Hv*t[1], mr*(Hu*n[0] + Hv*n[1])
            ti, ni = comp(R0*(1-1e-9)); to, no = comp(R0*(1+1e-9))
            okI &= np.allclose(ti, to, atol=1e-6) and np.allclose(ni, no, atol=1e-6)
    Hu, Hv, _, _ = f(np.array([400.0]), np.array([300.0]), 1000.0, 1.0, 2.0)
    okI &= abs(Hu[0]-1) < 1e-4 and abs(Hv[0]) < 1e-4
check("plotted fields: tangential H and normal B continuous at r=a,b; H -> H0 far away", okI)

# 15. shielding.py (one general-d implementation) equals an independent,
#     geometry-specific implementation written from the classical formulas
import shielding as sh
def ref_sphere(Z, X, mu, a, b):
    r = np.hypot(Z, X); r = np.where(r < 1e-12, 1e-12, r); c, s = Z/r, X/r
    k = (a/b)**3; den = (mu+2)*(2*mu+1) - 2*k*(mu-1)**2
    Hint = 9*mu/den; A = 3*(2*mu+1)/den; B = 3*(mu-1)*a**3/den
    C = b**3*(mu-1)*(2*mu+1)*(1-k)/den
    Hr = np.where(r < a, Hint*c, np.where(r <= b, (A-2*B/r**3)*c, (1+2*C/r**3)*c))
    Ht = np.where(r < a, -Hint*s, np.where(r <= b, -(A+B/r**3)*s, -(1-C/r**3)*s))
    return Hr*c - Ht*s, Hr*s + Ht*c, np.where((r >= a) & (r <= b), mu, 1.0), 1/Hint
def ref_cylinder(X, Y, mu, a, b):
    r = np.hypot(X, Y); r = np.where(r < 1e-12, 1e-12, r); c, s = X/r, Y/r
    D = (mu+1)**2*b**2 - (mu-1)**2*a**2; Hint = 4*mu*b**2/D
    E1 = -2*(mu+1)*b**2/D; F1 = -2*(mu-1)*a**2*b**2/D; D1 = (mu**2-1)*(b**2-a**2)*b**2/D
    Hr = np.where(r < a, Hint*c, np.where(r <= b, -(E1-F1/r**2)*c, (1+D1/r**2)*c))
    Ht = np.where(r < a, -Hint*s, np.where(r <= b, (E1+F1/r**2)*s, -(1-D1/r**2)*s))
    return Hr*c - Ht*s, Hr*s + Ht*c, np.where((r >= a) & (r <= b), mu, 1.0), 1/Hint
gU, gV = np.meshgrid(np.linspace(-4, 4, 121), np.linspace(-4, 4, 121))
ok15 = True
for mu_ in (0.2, 1.0, 7.0, 1000.0):
    for a_ in (0.3, 1.0, 1.7):
        for ref, d in ((ref_sphere, 3), (ref_cylinder, 2)):
            r1, r2 = ref(gU, gV, mu_, a_, 2.0), sh.fields(gU, gV, mu_, a_, 2.0, d)
            ok15 &= all(np.allclose(p, q, rtol=1e-11, atol=1e-13) for p, q in zip(r1, r2))
check("shielding.py fields = independent geometry-specific implementation", ok15)

# 16. the notebook carries exactly the core of shielding.py (text and numbers)
import json, re
nb = json.load(open("magnetic_shielding_interactive.ipynb"))
nb_core = next("".join(c["source"]) for c in nb["cells"]
               if c["cell_type"] == "code" and "def first_harmonic_coefficients" in "".join(c["source"]))
core = re.search(r"# ---- BEGIN CORE.*?\n(.*?)# ---- END CORE ----", open("shielding.py").read(), re.S).group(1)
ns = {"np": np}; exec(nb_core, ns)   # the notebook imports numpy in an earlier cell
same_num = all(np.allclose(p, q) for p, q in zip(ns["fields"](gU, gV, 1000.0, 1.0, 2.0, 3), sh.fields(gU, gV, 1000.0, 1.0, 2.0, 3)))
check("notebook core identical to shielding.py (text and results)", nb_core.strip() == core.strip() and same_num)

# 17. helper functions of shielding.py agree with the independent results above
ok17 = np.isclose(sh.X_STAR, (1 + np.sqrt(33))/16)
for m_ in (0.5, 10.0, 1000.0):
    for x in (0.2, 0.6):
        ok17 &= np.isclose(sh.SF_sph(m_, x), SF_sph(m_, x)) and np.isclose(sh.SF_cyl(m_, x), SF_cyl(m_, x))
        for n in (1, 2, 3):
            ok17 &= np.isclose(sh.harmonic_shielding_factor(m_, x, n, 2), S_cyl_nu(m_, x, 1.0, n))
            ok17 &= np.isclose(sh.harmonic_shielding_factor(m_, x, n, 3), S_sph_l(m_, x, 1.0, n))
        A_, B_, C_ = sh.three_term_coefficients(x, 3)
        ok17 &= np.isclose(A_*m_ + B_ + C_/m_, SF_sph(m_, x))
        ok17 &= np.isclose(sh.plus_one_approximation(m_, x, 3) - SF_sph(m_, x), A_*(2 - 1/m_))
for d, g in ((2, "cyl"), (3, "sph")):
    ok17 &= abs(cloak_D(sh.cloak_permeability(1.0, 1.3, d), 1.0, 1.3, g)) < 1e-12
ok17 &= np.isclose(sh.regime_crossover(0.5, 3), 9/(2*(1 - 0.125))) and np.isclose(sh.regime_crossover(0.5, 2), 4/(1 - 0.25))
check("shielding.py helper functions agree with independent solutions", ok17)


# ============================================================================
# 18-27. Statements of the final revision (main text, Supplementary Materials
#        I and II), each checked against the exact solution.
# ============================================================================
def SFd(mu, x, d):
    return 1 + (d - 1) / d**2 * (1 - x**d) * (mu - 1)**2 / mu

# 18. section 2.3: numbers quoted for the sphere with a/b = 1/2
A = (2/9) * (1 - 0.125)
mm = np.linspace(1, 10, 900001)
rel = (A*mm + 1)/SF_sph(mm, .5) - 1
check("2.3: '+1' excess 19% at mu=1, max 27% at mu~2.2, 14% at mu=10",
      round(100*rel[0]) == 19 and round(100*rel.max()) == 27
      and abs(mm[rel.argmax()] - 2.2) < 0.05 and round(100*rel[-1]) == 14)
hint = [100/SF_sph(m, .5) for m in (10, 100, 1000)]
check("2.3: cavity field 39%, 5.0%, 0.51% of H0 and SF ~ 2.6, 20, 195 (mu = 10, 100, 1000)",
      round(hint[0]) == 39 and round(hint[1], 1) == 5.0 and round(hint[2], 2) == 0.51
      and round(SF_sph(10, .5), 1) == 2.6 and round(SF_sph(100, .5)) == 20 and round(SF_sph(1000, .5)) == 195)
check("2.3: crossover mu* = 1/[C_d(1-k_d)] is about 5 for a/b = 1/2 (both shells)",
      round(9/(2*(1 - .125))) == 5 and round(4/(1 - .25)) == 5)

# 19. section 2.4: geometry comparison
xg = np.linspace(1e-4, 1 - 1e-4, 4001)
rat = (8/9)*(1 + xg + xg**2)/(1 + xg)
check("2.4: (SF_sph-1)/(SF_cyl-1) increases monotonically from 8/9 to 4/3",
      np.all(np.diff(rat) > 0) and abs(rat[0] - 8/9) < 1e-3 and abs(rat[-1] - 4/3) < 1e-3)
MU, XX = np.meshgrid(np.logspace(-4, 8, 121), np.linspace(0.001, 0.999, 120))
bracket = (2/9)*(1 - XX**3) - (1/4)*(1 - XX**2)
diff = SF_sph(MU, XX) - SF_cyl(MU, XX)
nz = np.abs(MU - 1) > 1e-9
check("2.4: SF_sph - SF_cyl has the sign of the geometric bracket for every mu != 1",
      np.all(np.sign(diff[nz]) == np.sign(bracket[nz])))
check("2.4: advantage of the cylinder stays below 9/8 and approaches it for x->0, mu->inf",
      np.all(SF_cyl(MU, XX)/SF_sph(MU, XX) < 9/8)
      and abs(SF_cyl(1e9, 1e-6)/SF_sph(1e9, 1e-6) - 9/8) < 1e-4)
check("2.4: advantage of the sphere tends to 4/3 for thin shells with mu(1-x) >> 1",
      abs(SF_sph(1e9, 0.999)/SF_cyl(1e9, 0.999) - 4/3) < 2e-3
      and abs(SF_sph(10, 0.999)/SF_cyl(10, 0.999) - 1) < 2e-3)
check("2.4: SF_sph/SF_cyl = 1.04 (a/b=0.5) and 0.95 (a/b=0.3, 216.8 vs 228.0) at mu=1000",
      round(SF_sph(1000, .5)/SF_cyl(1000, .5), 2) == 1.04 and round(SF_sph(1000, .3)/SF_cyl(1000, .3), 2) == 0.95
      and round(SF_sph(1000, .3), 1) == 216.8 and round(SF_cyl(1000, .3), 1) == 228.0)
check("2.4: C_d = N(1-N) with N = 1/d", all(np.isclose((d-1)/d**2, (1/d)*(1 - 1/d)) for d in (2, 3)))
# axial field on an infinite cylinder: the uniform field H0 z satisfies every interface condition
th_ = np.linspace(0, 2*np.pi, 37)
n_r = np.array([np.cos(th_), np.sin(th_), 0*th_])        # normal of the cylinder surfaces
H_ax = np.array([0, 0, 1.0])[:, None] * np.ones_like(th_)
check("2.4: a uniform axial field satisfies the interface conditions of an infinite cylinder (no shielding)",
      np.allclose((H_ax*n_r).sum(0), 0) and np.allclose(1000*(H_ax*n_r).sum(0), (H_ax*n_r).sum(0)))

# 20. section 2.5: the determinants reduce to the denominators of equations (1) and (2)
ok20 = True
for m_ in (0.3, 7.0, 1000.0):
    for x in (0.2, 0.7):
        ok20 &= np.isclose(det_c(m_, x, 1), (m_ + 1)**2 - (m_ - 1)**2*x**2)
        ok20 &= np.isclose(det_s(m_, x, 1), (m_ + 2)*(2*m_ + 1) - 2*x**3*(m_ - 1)**2)
check("2.5: nu = l = 1 determinants are the denominators of equations (1) and (2) (times b^2, b^3)", ok20)

# 21. section 3 and figure 2: fields for mu_r = 1000, a/b = 1/2
def field_at(rr, tt, mu, a, b, d):
    Hu, Hv, mr, SF = sh.fields(rr*np.cos(tt), rr*np.sin(tt), mu, a, b, d)
    return np.hypot(Hu, Hv), mr, SF
ok21a = ok21b = ok21c = ok21d = True
vals = {}
for d in (3, 2):
    a_, b_, mu_ = 1.0, 2.0, 1000.0
    Hint_ = 1/SFd(mu_, .5, d)
    rr, tt = np.meshgrid(np.linspace(a_*(1 + 1e-12), b_*(1 - 1e-12), 401), np.linspace(0, np.pi, 721))
    Hs, mr, _ = field_at(rr, tt, mu_, a_, b_, d)
    ok21a &= bool(np.all(mr == mu_))                                  # every sample lies in the wall
    Bs = mr*Hs
    i = np.unravel_index(Bs.argmax(), Bs.shape)
    ok21a &= np.isclose(rr[i], a_) and np.isclose(tt[i], np.pi/2) and np.isclose(Bs.max(), mu_*Hint_, rtol=1e-9)
    ok21a &= round(Bs.max()) == 5
    eps = 1e-10
    H_in_eq, _, _ = field_at(b_*(1 - eps), np.pi/2, mu_, a_, b_, d)
    H_out_eq, _, _ = field_at(b_*(1 + eps), np.pi/2, mu_, a_, b_, d)
    H_in_po, _, _ = field_at(b_*(1 - eps), 0.0, mu_, a_, b_, d)
    H_out_po, _, _ = field_at(b_*(1 + eps), 0.0, mu_, a_, b_, d)
    ok21b &= np.isclose(H_out_eq, H_in_eq, rtol=1e-6) and H_out_eq < 1e-2            # weak, = |H| in wall
    ok21b &= np.isclose(H_out_po, mu_*H_in_po, rtol=1e-6) and H_out_po > 1.5          # strong, = |B| in wall
    ok21c &= Hs.max() < 6e-3 and round(1e3*Hs.max()) == 5 and Hs.max() < 1e-2      # 'about 5e-3', > 2 orders below H0
    ok21c &= np.isclose(H_out_po/H_in_po, mu_, rtol=1e-6)                           # factor mu_r at the poles
    vals[d] = (SFd(mu_, .5, d), Hint_, H_out_po)
ok21d = (round(vals[3][0], 1) == 195.1 and round(vals[2][0], 1) == 188.1
         and round(1e3*vals[3][1], 1) == 5.1 and round(1e3*vals[2][1], 1) == 5.3)
check("3: max |B| in the wall = mu_r B_int, about 5 mu0 H0, on the inner surface at the equator", ok21a)
check("3: outside field = |H| of the wall at the equator (weak) and = |B| of the wall at the poles (strong)", ok21b)
check("3: |H| in the wall at most ~5e-3 H0, and smaller than in the vacuum at the poles by exactly mu_r", ok21c)
check("3/fig. 2: SF = 195.1 and 188.1, H_int/H0 = 5.1e-3 and 5.3e-3", ok21d)
mu_line = np.logspace(0, 5, 2001)
r05 = SF_sph(mu_line, .5)/SF_cyl(mu_line, .5)
check("3/fig. 3(a): for a/b = 1/2 the two curves differ by less than 4% at every mu_r",
      np.all((r05 >= 1 - 1e-12) & (r05 < 1.04)))
check("3: SF -> infinity as mu_r -> 0 (ideal superconducting shell)", SF_sph(1e-9, .5) > 1e7 and SF_cyl(1e-9, .5) > 1e7)

# 22. Supplementary Material I, section 5: higher harmonics
ok22 = True
for d in (2, 3):
    for x in np.linspace(0.3, 0.99, 50):     # x >= 0.3 keeps x**(m+s) resolvable in double precision
        Sm = [sh.harmonic_shielding_factor(1000.0, x, m, d) for m in range(1, 9)]
        ok22 &= bool(np.all(np.diff(Sm) > 0))
        for m in range(1, 5):
            ok22 &= np.isclose(sh.harmonic_shielding_factor(7.0, x, m, d),
                               sh.harmonic_shielding_factor(1/7.0, x, m, d))
check("SM I s5: S_m increases with m for every a/b, and S_m(mu) = S_m(1/mu)", ok22)
ok22b = True
for d in (2, 3):
    for m in (1, 2, 3):
        s_ = m + d - 2; t = 1e-6; mu_ = 1e9
        ok22b &= np.isclose((sh.harmonic_shielding_factor(mu_, 1 - t, m, d) - 1)/(mu_*t), m*s_/(m + s_), rtol=1e-4)
check("SM I s5: thin shell, S_m - 1 ~ [ms/(m+s)] mu_r t/b", ok22b)

# 23. Supplementary Material I, section 6
check("SM I s6: thin shells, (SF_sph-1)/(SF_cyl-1) -> 4/3 and SF ratio -> 4/3 only when mu t/b >> 1",
      abs((SF_sph(10, 1-1e-6) - 1)/(SF_cyl(10, 1-1e-6) - 1) - 4/3) < 1e-5
      and abs(SF_sph(1e12, 1-1e-6)/SF_cyl(1e12, 1-1e-6) - 4/3) < 1e-3)
ok23 = True
for d in (2, 3):
    N = 1/d
    for m_ in (0.5, 3.0, 50.0, 1e4):
        hint_, P_, Q_, R_ = sh.first_harmonic_coefficients(m_, 1e-30, 1.0, d)
        ok23 &= np.isclose(P_, 1/(1 + N*(m_ - 1)))                      # solid body
        hint_, P_, Q_, R_ = sh.first_harmonic_coefficients(m_, 0.4, 1.0, d)
        ok23 &= np.isclose(hint_/P_, m_/(1 + (1 - N)*(m_ - 1)))         # cavity amplification
check("SM I s6: solid-body factor 1/[1+N(mu-1)] and cavity factor mu/[1+(1-N)(mu-1)]", ok23)

# 24. Supplementary Material I, section 7: exercises
ok24 = True
for d in (2, 3):
    hint_, P_, Q_, R_ = sh.first_harmonic_coefficients(1e-12, 0.6, 1.0, d)
    ok24 &= abs(hint_) < 1e-9 and np.isclose(R_, -1.0/(d - 1), rtol=1e-9)
check("SM I s7, exercise 1: mu_r -> 0 gives H_int -> 0 and R -> -b^d H0/(d-1)", ok24)
ok24b = True
for d in (2, 3):
    a_, b_ = 0.7, 1.0
    mu_cl = ((d - 1)*b_**d + a_**d)/((d - 1)*(b_**d - a_**d))
    P_ = 1.0; Q_ = -P_*a_**d/(d - 1)
    dphi_a = -P_ + (1 - d)*Q_*a_**(-d)                 # d/dr of (-P r + Q r^(1-d)) at r = a
    ok24b &= abs(dphi_a) < 1e-12
    # with R = 0 outside, the two conditions at r = b fix P and mu_r
    M = np.array([[-b_ + b_**(1 - d)*(-a_**d/(d - 1))]]); P_sol = -b_/M[0, 0]
    dphi_b = P_sol*(-1 + (1 - d)*(-a_**d/(d - 1))*b_**(-d))
    ok24b &= np.isclose(mu_cl*dphi_b, -1.0)
check("SM I s7, exercise 2: Q = -P a^d/(d-1) and cloak permeability [(d-1)b^d+a^d]/[(d-1)(b^d-a^d)]", ok24b)

# 25. Supplementary Material II, activities 1-4
check("SM II act. 2: SF = 2.6, 20.1, 195.1 (fig. 2) and '+1' formula 2.94 vs exact 2.57 at mu=10",
      round(SF_sph(10, .5), 1) == 2.6 and round(SF_sph(100, .5), 1) == 20.1 and round(SF_sph(1000, .5), 1) == 195.1
      and round(A*10 + 1, 2) == 2.94 and round(SF_sph(10, .5), 2) == 2.57)
check("SM II act. 3: 109.2 vs 90.8 (a/b = 0.8) and 216.8 vs 228.0 (a/b = 0.3) at mu = 1000",
      round(SF_sph(1000, .8), 1) == 109.2 and round(SF_cyl(1000, .8), 1) == 90.8
      and round(SF_sph(1000, .3), 1) == 216.8 and round(SF_cyl(1000, .3), 1) == 228.0)
def bisect(m_):
    lo, hi = 0.01, 0.99
    f = lambda x: SF_sph(m_, x) - SF_cyl(m_, x)
    for _ in range(200):
        mid = 0.5*(lo + hi)
        lo, hi = (mid, hi) if f(lo)*f(mid) > 0 else (lo, mid)
    return 0.5*(lo + hi)
check("SM II act. 3: bisection gives a/b = 0.42 for mu = 10, 1000 and 1e4",
      all(abs(bisect(m_) - xstar) < 1e-9 and round(bisect(m_), 2) == 0.42 for m_ in (10, 1e3, 1e4)))
gu = np.linspace(-8, 8, 1601)
UU, VV = np.meshgrid(gu, gu)
Hu, Hv, mr, _ = sh.fields(UU, VV, 1000.0, 1.0, 2.0, 3)
Hm = np.hypot(Hu, Hv)
iH = np.unravel_index(Hm.argmax(), Hm.shape)
rH, tH = np.hypot(UU[iH], VV[iH]), np.arctan2(abs(VV[iH]), abs(UU[iH]))
H_pole_out = field_at(2.0*(1 + 1e-10), 0.0, 1000.0, 1.0, 2.0, 3)[0]
ring = field_at(2.0*(1 + 1e-10), np.linspace(0, np.pi, 1801), 1000.0, 1.0, 2.0, 3)[0]
check("SM II act. 1 and 4: |H| largest just outside the poles (about 3 H0); |B| outside the wall smallest near the equator",
      rH < 2.02 and tH < 0.05 and round(H_pole_out) == 3 and abs(np.linspace(0, np.pi, 1801)[ring.argmin()] - np.pi/2) < 1e-3)

print("\nALL CHECKS PASSED" if ok_all else "\nSOME CHECKS FAILED")
