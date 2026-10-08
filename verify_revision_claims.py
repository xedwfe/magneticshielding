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

print("\nALL CHECKS PASSED" if ok_all else "\nSOME CHECKS FAILED")
