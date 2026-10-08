"""
paper_figures.py -- reproduces every figure of the revised manuscript
"Magnetic shielding as a motivational tool in teaching classical
electromagnetic theory" (EJP-110956) and of its Supplementary Materials I
(figures 4 and 5) and II (figures 2 to 4).

Offline: numpy + matplotlib only. All figures are written as vector PDF at
their final printed width (446 pt = 6.2 in), so that fonts are legible at
print size. Colour maps are embedded as high-resolution raster layers inside
the PDF; all lines, arrows and text remain vector.

The exact fields come from shielding.py (Supplementary Material I,
section 4). Run from the repository folder:  python paper_figures.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Arc, Circle, FancyArrowPatch, Wedge

plt.rcParams.update({
    "font.size": 10, "axes.labelsize": 10.5, "axes.titlesize": 10.5,
    "xtick.labelsize": 9.5, "ytick.labelsize": 9.5, "legend.fontsize": 9.5,
    "mathtext.fontset": "cm", "font.family": "serif", "pdf.fonttype": 42,
})
TW = 6.2                     # text width in inches (iopart, 12 pt)
CMAP = "viridis"
NORM = LogNorm(vmin=1e-3, vmax=6.0)


# ------------------------------------------------- physics (shielding.py)
# The exact solution lives in shielding.py, the single source of the physics
# shared by all scripts of the repository.
from shielding import fields_sphere, fields_cylinder, SF_sph, SF_cyl  # noqa: E402


# ------------------------------------------------------------ field map
def field_map(ax, geom, mu, a, b, L=4.0, quantity="B", annotate=True,
              density=1.1, labels=("", "")):
    n = 520
    g = np.linspace(-L, L, n)
    U, V = np.meshgrid(g, g)
    if geom == "sph":
        Hu, Hv, murel, SF = fields_sphere(U, V, mu, a, b)
    else:
        Hu, Hv, murel, SF = fields_cylinder(U, V, mu, a, b)
    mag = np.hypot(Hu, Hv) * (murel if quantity == "B" else 1.0)
    im = ax.pcolormesh(U, V, np.clip(mag, 1.01e-3, 5.99), norm=NORM,
                       cmap=CMAP, shading="auto", rasterized=True)
    ax.streamplot(g, g, Hu * murel, Hv * murel, density=density, color="w",
                  linewidth=0.6, arrowsize=0.6, zorder=2)
    for R in (a, b):
        ax.add_patch(Circle((0, 0), R, fill=False, ec="k", lw=0.9, zorder=3))
    ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_aspect("equal")
    ax.set_xlabel(labels[0], labelpad=1); ax.set_ylabel(labels[1], labelpad=1)
    if annotate:
        m, e = f"{1/SF:.1e}".split("e")
        ax.text(0.03, 0.97, f"$SF = {SF:.1f}$\n$H_{{\\rm int}}/H_0 = {m}\\times10^{{{int(e)}}}$",
                transform=ax.transAxes, va="top", ha="left", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.3",
                          alpha=0.9), zorder=4)
    return im, SF


# ------------------------------------------------- figure 1: schematics
def schematic(ax, geom, polar=False):
    """Shell geometry; polar=True adds the polar coordinates used in SM I."""
    a, b = 1.0, 1.75
    ax.add_patch(Circle((0, 0), b, fc="#cfe0f5", ec="k", lw=1.2))
    ax.add_patch(Circle((0, 0), a, fc="white", ec="k", lw=1.2))
    for yy in (-1.25, -0.42, 0.42, 1.25):     # applied field arrows on the left
        ax.add_patch(FancyArrowPatch((-3.3, yy), (-2.25, yy),
                     arrowstyle="-|>", mutation_scale=11, lw=1.1, color="k"))
    ax.text(-2.8, 1.65, r"$\vec H_0$", ha="center", fontsize=11)
    ax.text(-2.8, -2.05, r"($\vec B_0=\mu_0\vec H_0$)", ha="center", fontsize=9)
    ang_a, ang_b = np.deg2rad(125), np.deg2rad(-35)
    ax.add_patch(FancyArrowPatch((0, 0), (a * np.cos(ang_a), a * np.sin(ang_a)),
                 arrowstyle="-|>", mutation_scale=9, lw=1.0, color="#1f4e9c"))
    ax.add_patch(FancyArrowPatch((0, 0), (b * np.cos(ang_b), b * np.sin(ang_b)),
                 arrowstyle="-|>", mutation_scale=9, lw=1.0, color="#1f4e9c"))
    ax.text(-0.52, 0.33, "$a$", fontsize=11, color="#1f4e9c")
    ax.text(0.75, -0.42, "$b$", fontsize=11, color="#1f4e9c")
    ax.text(0.72, 1.12, r"$\mu_r$", ha="center", fontsize=11)
    ax.text(-0.42, -0.62, r"$\mu_0$", ha="center", fontsize=10)
    ax.text(2.25, 1.45, r"$\mu_0$", ha="center", fontsize=10)
    ax.annotate("", xy=(3.1, 0), xytext=(-3.4, 0),
                arrowprops=dict(arrowstyle="-|>", lw=0.7, color="0.35",
                                linestyle=(0, (4, 3))))
    if geom == "sph":
        ax.text(3.1, 0.15, "$z$", fontsize=11)
        ax.text(0, -2.3, "axisymmetric about $z$", ha="center", fontsize=9)
    else:
        ax.text(3.1, 0.15, "$x$", fontsize=11)
        ax.annotate("", xy=(0, 2.5), xytext=(0, -2.2),
                    arrowprops=dict(arrowstyle="-|>", lw=0.7, color="0.35",
                                    linestyle=(0, (4, 3))))
        ax.text(0.15, 2.4, "$y$", fontsize=11)
        ax.add_patch(Circle((0, 0), 0.09, fc="white", ec="k", lw=0.8, zorder=5))
        ax.add_patch(Circle((0, 0), 0.025, fc="k", ec="k", zorder=6))
        ax.text(0.14, 0.10, "$z$", fontsize=10)
        ax.text(0, -2.3, "infinitely long along $z$", ha="center", fontsize=9)
    if polar:                                  # (rho, theta) or (r, theta), theta from H0
        ang, R = np.deg2rad(20), 2.45
        ax.add_patch(FancyArrowPatch((0, 0), (R * np.cos(ang), R * np.sin(ang)),
                     arrowstyle="-|>", mutation_scale=9, lw=1.0, color="#2a8c3a", zorder=4))
        ax.text(R * np.cos(ang) + 0.02, R * np.sin(ang) + 0.12,
                r"$\rho$" if geom == "cyl" else "$r$", fontsize=11, color="#2a8c3a")
        ax.add_patch(Arc((0, 0), 1.3, 1.3, theta1=0, theta2=20, color="0.25", lw=0.9, zorder=4))
        ax.text(0.74, 0.07, r"$\theta$", fontsize=10)
    ax.set_xlim(-3.5, 3.4); ax.set_ylim(-2.6, 2.6); ax.set_aspect("equal")
    ax.axis("off")


def figure1():
    fig, axs = plt.subplots(1, 2, figsize=(TW, 2.55))
    schematic(axs[0], "sph"); schematic(axs[1], "cyl")
    axs[0].set_title("(a) spherical shell (section through $z$)", pad=2)
    axs[1].set_title("(b) cylindrical shell (cross-section)", pad=2)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.9, bottom=0.01, wspace=0.05)
    fig.savefig("fig1_geometry.pdf")
    for i, g in enumerate(("sph", "cyl")):     # single panels for SM I
        f, ax = plt.subplots(figsize=(3.4, 2.6)); schematic(ax, g, polar=True)
        f.subplots_adjust(left=0, right=1, top=1, bottom=0)
        f.savefig(f"sm1_{'sphere' if g == 'sph' else 'cylinder'}_section.pdf")
        plt.close(f)
    plt.close(fig)


# ------------------------------------------------ figure 2: field maps
def figure2(mu=1000.0, a=1.0, b=2.0):
    fig, axs = plt.subplots(1, 2, figsize=(TW, 3.05), constrained_layout=True)
    im, _ = field_map(axs[0], "sph", mu, a, b, labels=("$z$ (m)", "$x$ (m)"))
    field_map(axs[1], "cyl", mu, a, b, labels=("$x$ (m)", "$y$ (m)"))
    axs[0].set_title(r"(a) sphere, meridian plane", pad=3)
    axs[1].set_title(r"(b) cylinder, transverse plane", pad=3)
    cb = fig.colorbar(im, ax=axs, shrink=0.92, pad=0.01)
    cb.set_label(r"$|\vec B|/(\mu_0 H_0)$")
    fig.savefig("fig2_field_maps.pdf", dpi=300)
    plt.close(fig)


# ------------------------------------- figure 3: SF(mu) and the geometry
def figure3(a=1.0, b=2.0):
    x0 = a / b
    fig = plt.figure(figsize=(TW, 3.15))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1], height_ratios=[4.2, 1],
                          hspace=0.06, wspace=0.32, left=0.095, right=0.985,
                          top=0.92, bottom=0.14)
    ax = fig.add_subplot(gs[0, 0]); sx = fig.add_subplot(gs[1, 0], sharex=ax)
    bx = fig.add_subplot(gs[:, 1])
    mu = np.logspace(0, 5, 600)
    ax.loglog(mu, SF_sph(mu, x0), "-", color="#1f4e9c", lw=1.8, label="sphere")
    ax.loglog(mu, SF_cyl(mu, x0), "--", color="#2a8c3a", lw=1.8, label="cylinder")
    ax.loglog(mu, (2/9)*(1-x0**3)*mu, ":", color="0.45", lw=1.0)
    ax.axhline(1, color="0.6", lw=0.7)
    ax.annotate("plateau, $SF\\approx 1$", xy=(1.7, 1.08), xytext=(1.15, 40),
                fontsize=8.8, arrowprops=dict(arrowstyle="-", lw=0.6))
    ax.annotate(r"$SF\simeq \frac{2}{9}(1-k_3)\,\mu_r$", xy=(4e3, 8.6e2),
                xytext=(5e2, 12), fontsize=9,
                arrowprops=dict(arrowstyle="-", lw=0.6))
    ax.set_xlim(1, 1e5); ax.set_ylim(0.8, 3e4)
    ax.set_ylabel(r"$SF = H_0/H_{\rm int}$")
    ax.set_title(f"(a) $a/b = {x0:.1f}$", pad=3)
    ax.legend(loc="upper left", frameon=False)
    ax.grid(True, which="major", alpha=0.25)
    plt.setp(ax.get_xticklabels(), visible=False)
    for i, (lo, hi, lab) in enumerate([(1e2, 1e4, "soft ferrites"),
                                       (1e3, 1e4, "silicon steels"),
                                       (1e4, 1e5, "Ni-Fe alloys")]):
        y = 2 - i
        sx.plot([lo, hi], [y, y], color="#b4522c", lw=3.0, solid_capstyle="butt")
        sx.text(lo / 1.3, y, lab, ha="right", va="center", fontsize=8.4,
                color="#7a3317")
    sx.set_ylim(-0.7, 2.7); sx.set_yticks([])
    sx.set_xscale("log"); sx.set_xlim(1, 1e5)
    sx.set_xlabel(r"relative permeability $\mu_r$")
    for sp in ("top", "right", "left"):
        sx.spines[sp].set_visible(False)

    x = np.linspace(0.0, 0.999, 800)
    xs = (1 + np.sqrt(33)) / 16
    bx.axvspan(0, xs, color="#2a8c3a", alpha=0.07, lw=0)
    bx.axvspan(xs, 1, color="#1f4e9c", alpha=0.07, lw=0)
    for m, ls in [(10, ":"), (100, "-."), (1e3, "--"), (1e5, "-")]:
        bx.plot(x, SF_sph(m, x) / SF_cyl(m, x), ls, lw=1.3, color="k",
                alpha=0.35 + 0.13 * np.log10(m),
                label=rf"$\mu_r=10^{{{int(np.log10(m))}}}$")
    env = (8 / 9) * (1 + x + x ** 2) / (1 + x)
    bx.plot(x, env, color="#b4522c", lw=1.2, label=r"$\mu_r\to\infty$")
    bx.axhline(1, color="0.6", lw=0.7)
    bx.axvline(xs, color="0.3", lw=0.8, ls=(0, (2, 2)))
    bx.text(xs + 0.015, 0.855, r"$x^\ast\approx0.42$", fontsize=9)
    bx.text(0.03, 1.37, "cylinder\nbetter", fontsize=8.8, color="#1d6b2b", va="top")
    bx.text(0.46, 1.37, "sphere better", fontsize=8.8, color="#1f4e9c", va="top")
    bx.annotate("4/3", xy=(0.997, 4 / 3), xytext=(0.70, 1.29), fontsize=9,
                arrowprops=dict(arrowstyle="-", lw=0.6))
    bx.annotate("8/9", xy=(0.0, 8 / 9), xytext=(0.08, 0.86), fontsize=9,
                arrowprops=dict(arrowstyle="-", lw=0.6))
    bx.set_xlim(0, 1); bx.set_ylim(0.84, 1.38)
    bx.set_xlabel(r"$x = a/b$"); bx.set_ylabel(r"$SF_{\rm sph}/SF_{\rm cyl}$")
    bx.set_title("(b) sphere versus cylinder", pad=3)
    bx.legend(loc="lower right", frameon=False, fontsize=8.6)
    fig.savefig("fig3_shielding_factor.pdf")
    plt.close(fig)


# --------------------------------------------- Supplementary Material II
def sm2_mu_sweep(a=1.0, b=2.0):
    fig, axs = plt.subplots(1, 3, figsize=(TW, 2.45), constrained_layout=True)
    for ax, mu in zip(axs, (10, 100, 1000)):
        im, SF = field_map(ax, "sph", mu, a, b, annotate=False, density=0.8,
                           labels=("$z$ (m)", "$x$ (m)" if mu == 10 else ""))
        ax.set_title(rf"$\mu_r = {mu}$:  $SF = {SF:.1f}$", pad=3)
    cb = fig.colorbar(im, ax=axs, shrink=0.95, pad=0.01)
    cb.set_label(r"$|\vec B|/(\mu_0 H_0)$")
    fig.savefig("sm2_fig2_mu_sweep.pdf", dpi=300)
    plt.close(fig)


def sm2_geometry(mu=1000.0, b=2.0):
    fig, axs = plt.subplots(2, 2, figsize=(TW, 5.1), constrained_layout=True)
    for i, x in enumerate((0.8, 0.3)):
        for j, geom in enumerate(("sph", "cyl")):
            ax = axs[i, j]
            im, SF = field_map(ax, geom, mu, x * b, b, annotate=False,
                               density=0.8, labels=(
                                   ("$z$ (m)" if geom == "sph" else "$x$ (m)") if i else "",
                                   ("$x$ (m)" if geom == "sph" else "$y$ (m)")))
            name = "sphere" if geom == "sph" else "cylinder"
            ax.set_title(rf"{name}, $a/b = {x}$:  $SF = {SF:.1f}$", pad=3)
    cb = fig.colorbar(im, ax=axs, shrink=0.6, pad=0.01)
    cb.set_label(r"$|\vec B|/(\mu_0 H_0)$")
    fig.savefig("sm2_fig3_geometry.pdf", dpi=300)
    plt.close(fig)


def sm2_B_vs_H(mu=1000.0, a=1.0, b=2.0):
    fig, axs = plt.subplots(1, 2, figsize=(TW, 3.05), constrained_layout=True)
    im, _ = field_map(axs[0], "sph", mu, a, b, quantity="B", annotate=False,
                      labels=("$z$ (m)", "$x$ (m)"))
    field_map(axs[1], "sph", mu, a, b, quantity="H", annotate=False,
              labels=("$z$ (m)", ""))
    axs[0].set_title(r"(a) $|\vec B|/(\mu_0H_0)$", pad=3)
    axs[1].set_title(r"(b) $|\vec H|/H_0$", pad=3)
    cb = fig.colorbar(im, ax=axs, shrink=0.92, pad=0.01)
    cb.set_label("field magnitude (normalized)")
    fig.savefig("sm2_fig4_B_vs_H.pdf", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    figure1(); figure2(); figure3()
    sm2_mu_sweep(); sm2_geometry(); sm2_B_vs_H()
    # numbers quoted in the text
    Z, X = np.meshgrid(np.linspace(-4, 4, 801), np.linspace(-4, 4, 801))
    Hz, Hx, mr, SF = fields_sphere(Z, X, 1000.0, 1.0, 2.0)
    shell = mr > 1
    Bmax = (np.hypot(Hz, Hx) * mr)[shell].max()
    Hshell = np.hypot(Hz, Hx)[shell]
    print(f"sphere mu=1000, a/b=1/2: SF={SF:.2f}; max |B|/mu0H0 in shell = {Bmax:.2f}; "
          f"|H|/H0 in shell in [{Hshell.min():.1e}, {Hshell.max():.1e}]")
    print("figures written: fig1_geometry.pdf fig2_field_maps.pdf fig3_shielding_factor.pdf "
          "sm1_sphere_section.pdf sm1_cylinder_section.pdf sm2_fig2_mu_sweep.pdf "
          "sm2_fig3_geometry.pdf sm2_fig4_B_vs_H.pdf")
