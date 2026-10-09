"""
build_notebook.py: writes magnetic_shielding_interactive.ipynb.

The physics block of the notebook is copied verbatim from shielding.py
(between the CORE markers), so that the notebook is self-contained for
Google Colab and yet identical to the module used by the other scripts.
Run this script again whenever shielding.py changes.
"""
import json
import re

src = open("shielding.py").read()
core = re.search(r"# ---- BEGIN CORE.*?\n(.*?)# ---- END CORE ----", src, re.S).group(1)

md = lambda text: {"cell_type": "markdown", "metadata": {}, "source": text.strip("\n")}
code = lambda text: {"cell_type": "code", "metadata": {}, "execution_count": None,
                     "outputs": [], "source": text.strip("\n")}

cells = [
md(r"""
# Magnetic shielding by permeable shells: interactive notebook

Companion notebook to the article *Magnetic shielding as a context for teaching
magnetostatics in matter: spherical and cylindrical shells with open code*.

It plots the **exact** static field of a permeable spherical shell and of an
infinitely long cylindrical shell in a uniform applied field $H_0$, and lets you
change the parameters with sliders. It runs in a web browser on Google Colab or
Binder; nothing needs to be installed.

**How to use.** Run all cells (on Colab: *Runtime > Run all*; on Jupyter or
Binder: *Run > Run All Cells*) and move the sliders in section 2.

**What the maps show.** The colour gives $|\vec B|/(\mu_0H_0)$ on a logarithmic
scale. The white curves are streamlines of $\vec B$, everywhere tangent to the
field; their spacing carries no information about the field strength, which is
read from the colour alone. The boxes give the cavity field $H_{\rm int}/H_0$ and
the shielding factor $SF = H_0/H_{\rm int}$.

**Model.** Linear, isotropic and homogeneous shell of relative permeability
$\mu_r$, inner radius $a$, outer radius $b$, in vacuum; $\vec B = \mu_0\vec H$ in
vacuum and $\vec B = \mu_0\mu_r\vec H$ in the shell, whose wall is the region
$a<r<b$. The sphere is shown in a meridional plane and the cylinder in its
transverse plane. In both panels the applied field points to the right, along
the axis labelled $x$; for the sphere this axis is the symmetry axis, called $z$
in the article and in its supplements. The angle $\theta$ is measured from the
direction of the applied field: the poles of each surface are its points with
$\theta = 0$ and $\theta = \pi$, and its equator is formed by the points with
$\theta = \pi/2$.
"""),
code(r"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle
import ipywidgets as widgets
from ipywidgets import interact
%matplotlib inline
"""),
md(r"""
## 1. The exact solution

Both shells are solved at once, with the dimension $d$ of the problem as a
parameter ($d=3$ for the sphere, $d=2$ for the cylinder; Supplementary
Material I). The shielding factor is
$$SF = 1 + C_d\,(1-k_d)\,\frac{(\mu_r-1)^2}{\mu_r},\qquad C_d = \frac{d-1}{d^2},\qquad k_d = \left(\frac{a}{b}\right)^d.$$
The next cell is identical to the core of `shielding.py` in the repository.
"""),
code(core),
md(r"""
## 2. Explore the two shells (Activities 1 to 3)

Sliders: relative permeability `mu_r`, radius ratio `a/b`, applied field `H0`,
outer radius `b`, half-width `L` of the window and grid resolution `N`.
Changing `H0` leaves every normalized quantity unchanged: the model is linear, so
the shielding factor does not depend on the applied field (in real materials it
does, through saturation).
"""),
code(r"""
NORM = LogNorm(vmin=1e-3, vmax=6.0)

def sci(v):
    # number for the text boxes: plain if moderate, otherwise m x 10^e
    e = int(np.floor(np.log10(abs(v))))
    if -2 < e < 3:
        return f"{v:.3g}"
    return f"{v / 10**e:.2f}\\times10^{{{e}}}"

def draw_maps(mu_r=1000.0, ratio=0.5, H0=1.0, b=2.0, L=4.0, N=350, quantity="B"):
    a = ratio * b
    g = np.linspace(-L, L, int(N))
    U, V = np.meshgrid(g, g)
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.4))
    panels = ((axs[0], 3, "Spherical shell (meridional plane)"),
              (axs[1], 2, "Cylindrical shell (transverse plane)"))
    for ax, d, name in panels:
        Hu, Hv, mu_map, SF = fields(U, V, mu_r, a, b, d, H0)
        weight = mu_map if quantity == "B" else 1.0
        mag = np.hypot(Hu, Hv) * weight / H0
        im = ax.pcolormesh(U, V, np.clip(mag, 1.01e-3, 5.99), norm=NORM,
                           cmap="viridis", shading="auto")
        ax.streamplot(g, g, Hu * mu_map, Hv * mu_map, density=1.2, color="w",
                      linewidth=0.7, arrowsize=0.8)
        for R in (a, b):
            ax.add_patch(Circle((0, 0), R, fill=False, ec="k", lw=1.2))
        ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_aspect("equal")
        ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)")
        ax.set_title(f"{name},  $\\mu_r = {mu_r:.4g}$")
        ax.text(0.02, 0.97, f"$H_{{\\rm int}}/H_0 = {sci(1/SF)}$\n$SF = {SF:.1f}$",
                transform=ax.transAxes, va="top",
                bbox=dict(boxstyle="round", fc="white", alpha=0.9))
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label(r"$|\vec B|/(\mu_0 H_0)$ (log)" if quantity == "B"
                     else r"$|\vec H|/H_0$ (log)")
    plt.tight_layout()
    plt.show()

interact(draw_maps,
         mu_r=widgets.FloatLogSlider(value=1000, base=10, min=0, max=5, step=0.01, description="mu_r"),
         ratio=widgets.FloatSlider(value=0.5, min=0.05, max=0.95, step=0.01, description="a/b"),
         H0=widgets.FloatSlider(value=1.0, min=0.1, max=5.0, step=0.1, description="H0"),
         b=widgets.FloatSlider(value=2.0, min=0.5, max=3.5, step=0.1, description="b"),
         L=widgets.FloatSlider(value=4.0, min=2.0, max=8.0, step=0.5, description="L"),
         N=widgets.IntSlider(value=350, min=100, max=600, step=10, description="N"),
         quantity=widgets.fixed("B"));
"""),
md(r"""
## 3. The shielding factor as a function of $\mu_r$ (Activity 2)

The exact curves of article figure 3(a), for the value of $a/b$ chosen with the
slider. The dotted line is the common high-permeability approximation
$SF \approx C_d(1-k_d)\mu_r + 1$ for the sphere, which exceeds the exact value by
$C_d(1-k_d)(2-1/\mu_r)$ (in article figure 3(a) the dotted line is instead the
linear behaviour $C_d(1-k_d)\mu_r$, which meets the plateau at the crossover).
"""),
code(r"""
def plot_sf(ratio=0.5):
    mu = np.logspace(0, 5, 500)
    fig, ax = plt.subplots(figsize=(7, 4.8))
    ax.loglog(mu, shielding_factor(mu, ratio, 3), "-", lw=2, label="sphere")
    ax.loglog(mu, shielding_factor(mu, ratio, 2), "--", lw=2, label="cylinder")
    A = 2 / 9 * (1 - ratio ** 3)
    ax.loglog(mu, A * mu + 1, ":", color="0.4", label="sphere, '+1' approximation")
    ax.axhline(1, color="0.6", lw=0.8)
    ax.set_xlabel(r"relative permeability $\mu_r$"); ax.set_ylabel(r"$SF = H_0/H_{\rm int}$")
    ax.set_title(f"exact shielding factor, a/b = {ratio:.2f}")
    ax.legend(); ax.grid(True, which="major", alpha=0.3)
    plt.show()
    for m in (10, 100, 1000):
        s = shielding_factor(m, ratio, 3)
        print(f"mu_r = {m:>5}:  SF_sph = {s:8.2f}   H_int/H0 = {100/s:6.2f} %"
              f"   SF_cyl = {shielding_factor(m, ratio, 2):8.2f}")

interact(plot_sf, ratio=widgets.FloatSlider(value=0.5, min=0.05, max=0.95, step=0.01, description="a/b"));
"""),
md(r"""
## 4. Sphere or cylinder? (Activity 3)

Choose $\mu_r$ and $a/b$ and compare the two shells. Then find, by moving the
`a/b` slider, the value at which the two shielding factors are equal, and repeat
for other permeabilities. The cell after this one checks your answer.
"""),
code(r"""
def compare(mu_r=1000.0, ratio=0.3):
    s, c = shielding_factor(mu_r, ratio, 3), shielding_factor(mu_r, ratio, 2)
    better = "sphere" if s > c else ("cylinder" if c > s else "neither")
    print(f"SF_sph = {s:.2f}   SF_cyl = {c:.2f}   better shield: {better}")

interact(compare,
         mu_r=widgets.FloatLogSlider(value=1000, base=10, min=0, max=5, step=0.01, description="mu_r"),
         ratio=widgets.FloatSlider(value=0.3, min=0.05, max=0.95, step=0.005, description="a/b"));
"""),
code(r"""
# Check your answer: crossing point by bisection, compared with the exact value
def crossing_point(mu_r, lo=0.01, hi=0.99, tol=1e-12):
    f = lambda x: shielding_factor(mu_r, x, 3) - shielding_factor(mu_r, x, 2)
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(lo) * f(mid) > 0 else (lo, mid)
    return 0.5 * (lo + hi)

for m in (10, 1000, 1e4):
    print(f"mu_r = {m:>7g}:  a/b = {crossing_point(m):.6f}")
print(f"exact: (1 + sqrt(33))/16 = {X_STAR:.6f}, independent of mu_r")
"""),
md(r"""
## 5. $\vec B$ versus $\vec H$ (Activity 4)

The same configuration shown as $|\vec B|/(\mu_0H_0)$ and as $|\vec H|/H_0$, on
the same colour scale. Where is each one largest, and why? Why is the field just
outside the shell weak near the equator and strong near the poles? (Hint: which
component of which field is continuous at the surface?)
"""),
code(r"""
def draw_B_vs_H(mu_r=1000.0, ratio=0.5, geometry="sphere"):
    d = 3 if geometry == "sphere" else 2
    b, L = 2.0, 4.0
    g = np.linspace(-L, L, 350)
    U, V = np.meshgrid(g, g)
    Hu, Hv, mu_map, SF = fields(U, V, mu_r, ratio * b, b, d)
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.2))
    for ax, weight, label in ((axs[0], mu_map, r"$|\vec B|/(\mu_0H_0)$"),
                              (axs[1], 1.0, r"$|\vec H|/H_0$")):
        im = ax.pcolormesh(U, V, np.clip(np.hypot(Hu, Hv) * weight, 1.01e-3, 5.99),
                           norm=NORM, cmap="viridis", shading="auto")
        ax.streamplot(g, g, Hu * mu_map, Hv * mu_map, density=1.1, color="w",
                      linewidth=0.6, arrowsize=0.7)
        for R in (ratio * b, b):
            ax.add_patch(Circle((0, 0), R, fill=False, ec="k", lw=1.2))
        ax.set_aspect("equal"); ax.set_title(f"{geometry}: {label}")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    plt.tight_layout(); plt.show()
    shell = mu_map > 1
    print(f"in the shell: max |B|/(mu0 H0) = {(np.hypot(Hu, Hv) * mu_map)[shell].max():.2f},"
          f"  max |H|/H0 = {np.hypot(Hu, Hv)[shell].max():.2e}")

interact(draw_B_vs_H,
         mu_r=widgets.FloatLogSlider(value=1000, base=10, min=0, max=5, step=0.01, description="mu_r"),
         ratio=widgets.FloatSlider(value=0.5, min=0.05, max=0.95, step=0.01, description="a/b"),
         geometry=widgets.Dropdown(options=["sphere", "cylinder"], value="sphere"));
"""),
md(r"""
## 6. Exercises for students who wish to modify the code

1. Plot $H_{\rm int}/H_0$ instead of $SF$ as a function of $\mu_r$.
2. Allow $\mu_r < 1$ (change `min` of the `mu_r` slider to a negative exponent)
   and check that the shielding factor is unchanged under $\mu_r \to 1/\mu_r$.
   Interpret the limit $\mu_r \to 0$ (Supplementary Material I, section 6).
3. Verify numerically the cloaking condition of Supplementary Material I,
   section 7 (exercise 2): for a superconducting layer of outer radius $a$
   lining a shell of outer radius $b$ and relative permeability
   $\mu_r = [(d-1)b^d + a^d]/[(d-1)(b^d - a^d)]$, the field outside the shell
   equals the applied field.
4. Add a second concentric shell; this requires solving a larger linear system
   with the same interface conditions.
"""),
]

nb = {"cells": cells,
      "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                   "language_info": {"name": "python"}},
      "nbformat": 4, "nbformat_minor": 5}
for i, c in enumerate(nb["cells"]):
    c["id"] = f"cell-{i:02d}"
    c["source"] = [line for line in c["source"].splitlines(keepends=True)]
json.dump(nb, open("magnetic_shielding_interactive.ipynb", "w"), indent=1)
print("written magnetic_shielding_interactive.ipynb with", len(cells), "cells")
