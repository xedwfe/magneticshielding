"""
explore_shielding.py: explore the exact solution without Jupyter.

Replaces magnetic_shieldingV3_corrected.py and
magnetic_shielding_S_graf_corrected.py. Every command prints the relevant
numbers; the plotting commands also save a PNG file (choose its name with
--out, add --show to open a window).

Examples (one per activity of Supplementary Material II)
  python explore_shielding.py maps --mu 1000 --ratio 0.5           # Activity 1
  python explore_shielding.py table --ratio 0.5 --mus 1 2 10 100 1000  # Activity 2
  python explore_shielding.py curve --ratio 0.5                    # Activity 2 / SF curve
  python explore_shielding.py compare --mu 1000 --ratios 0.8 0.3   # Activity 3
  python explore_shielding.py crossover --mus 10 1000 10000        # Activity 3
  python explore_shielding.py maps --mu 1000 --ratio 0.5 --quantity H   # Activity 4

Offline: numpy and matplotlib only.
"""
import argparse
import numpy as np
import matplotlib

from shielding import (fields, shielding_factor, SF_sph, SF_cyl, X_STAR,
                       plus_one_approximation, regime_crossover,
                       ratio_envelope)

NAMES = {3: "sphere", 2: "cylinder"}


# ------------------------------------------------------------- plotting
def _plt(show):
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _map(ax, d, mu_r, a, b, L, n, quantity):
    from matplotlib.colors import LogNorm
    from matplotlib.patches import Circle
    g = np.linspace(-L, L, n)
    U, V = np.meshgrid(g, g)
    Hu, Hv, mu_map, SF = fields(U, V, mu_r, a, b, d)
    mag = np.hypot(Hu, Hv) * (mu_map if quantity == "B" else 1.0)
    im = ax.pcolormesh(U, V, np.clip(mag, 1.01e-3, 5.99),
                       norm=LogNorm(vmin=1e-3, vmax=6.0), cmap="viridis",
                       shading="auto")
    ax.streamplot(g, g, Hu * mu_map, Hv * mu_map, density=1.1, color="w",
                  linewidth=0.6, arrowsize=0.7)
    for R in (a, b):
        ax.add_patch(Circle((0, 0), R, fill=False, ec="k", lw=1.0))
    ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_aspect("equal")
    ax.set_xlabel("z (m)" if d == 3 else "x (m)")
    ax.set_ylabel("x (m)" if d == 3 else "y (m)")
    ax.text(0.02, 0.97, f"SF = {SF:.1f}\nH_int/H0 = {1/SF:.2e}",
            transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(boxstyle="round", fc="white", alpha=0.9))
    return im, SF


def cmd_maps(args):
    plt = _plt(args.show)
    a, b = args.ratio * args.b, args.b
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.8), constrained_layout=True)
    for ax, d in zip(axs, (3, 2)):
        im, SF = _map(ax, d, args.mu, a, b, args.L, args.N, args.quantity)
        ax.set_title(f"{NAMES[d]}, mu_r = {args.mu:g}, a/b = {args.ratio:g}")
        print(f"{NAMES[d]:8s}: SF = {SF:9.3f}   H_int/H0 = {1/SF:.3e}")
    label = "|B|/(mu0 H0)" if args.quantity == "B" else "|H|/H0"
    fig.colorbar(im, ax=axs, shrink=0.9, label=label + " (log)")
    _finish(plt, fig, args, f"maps_{args.quantity}.png")


def cmd_curve(args):
    plt = _plt(args.show)
    x = args.ratio
    mu = np.logspace(0, 5, 500)
    fig, ax = plt.subplots(figsize=(6.5, 4.6), constrained_layout=True)
    for d, ls in ((3, "-"), (2, "--")):
        ax.loglog(mu, shielding_factor(mu, x, d), ls, lw=2, label=NAMES[d])
        print(f"{NAMES[d]:8s}: plateau-to-linear crossover near mu_r = "
              f"{regime_crossover(x, d):.2f}")
    ax.axhline(1, color="0.6", lw=0.8)
    ax.set_xlabel("relative permeability mu_r")
    ax.set_ylabel("SF = H0/H_int")
    ax.set_title(f"exact shielding factor, a/b = {x:g}")
    ax.legend(); ax.grid(True, which="major", alpha=0.3)
    _finish(plt, fig, args, "shielding_factor.png")


def cmd_compare(args):
    plt = _plt(args.show)
    n = len(args.ratios)
    fig, axs = plt.subplots(n, 2, figsize=(9, 4.3 * n), squeeze=False,
                            constrained_layout=True)
    for i, x in enumerate(args.ratios):
        s, c = SF_sph(args.mu, x), SF_cyl(args.mu, x)
        better = "sphere" if s > c else ("cylinder" if c > s else "equal")
        print(f"a/b = {x:4.2f}: SF_sph = {s:9.2f}  SF_cyl = {c:9.2f}  "
              f"better: {better}")
        for j, d in enumerate((3, 2)):
            im, SF = _map(axs[i, j], d, args.mu, x * args.b, args.b, args.L,
                          args.N, "B")
            axs[i, j].set_title(f"{NAMES[d]}, a/b = {x:g}: SF = {SF:.1f}")
    fig.colorbar(im, ax=axs, shrink=0.6, label="|B|/(mu0 H0) (log)")
    _finish(plt, fig, args, "compare.png")


def crossing_point(mu_r, lo=0.01, hi=0.99, tol=1e-12):
    """a/b at which SF_sph = SF_cyl, by bisection."""
    f = lambda x: SF_sph(mu_r, x) - SF_cyl(mu_r, x)
    if f(lo) * f(hi) > 0:
        return float("nan")
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(lo) * f(mid) > 0 else (lo, mid)
    return 0.5 * (lo + hi)


def cmd_crossover(args):
    print(f"exact value: x* = (1 + sqrt(33))/16 = {X_STAR:.10f}")
    for m in args.mus:
        print(f"mu_r = {m:>9g}: bisection gives a/b = {crossing_point(m):.10f}")
    print("for a/b < x* the cylinder shields better; for a/b > x* the sphere,"
          f" by at most 4/3 (ratio envelope at a/b = 0.99: "
          f"{ratio_envelope(0.99):.3f})")


def cmd_table(args):
    x = args.ratio
    print(f"a/b = {x:g}")
    print(f"{'mu_r':>9} {'SF_sph':>10} {'SF_cyl':>10} {'H_int/H0 (sph)':>15}"
          f" {'+1 approx (sph)':>16} {'excess':>8}")
    for m in args.mus:
        s, c = SF_sph(m, x), SF_cyl(m, x)
        ap = plus_one_approximation(m, x, 3)
        print(f"{m:>9g} {s:>10.3f} {c:>10.3f} {1/s:>15.4e} {ap:>16.3f}"
              f" {100*(ap/s - 1):>7.1f}%")


def _finish(plt, fig, args, default):
    out = args.out or default
    fig.savefig(out, dpi=150)
    print(f"figure saved to {out}")
    if args.show:
        plt.show()


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    def common(sp, plot=True):
        sp.add_argument("--b", type=float, default=2.0, help="outer radius (m)")
        if plot:
            sp.add_argument("--L", type=float, default=4.0, help="half-width of the window (m)")
            sp.add_argument("--N", type=int, default=350, help="grid points per axis")
            sp.add_argument("--out", default=None, help="output PNG file")
            sp.add_argument("--show", action="store_true", help="open a window")

    sp = sub.add_parser("maps", help="field maps of both shells")
    sp.add_argument("--mu", type=float, default=1000.0)
    sp.add_argument("--ratio", type=float, default=0.5, help="a/b")
    sp.add_argument("--quantity", choices=("B", "H"), default="B")
    common(sp); sp.set_defaults(func=cmd_maps)

    sp = sub.add_parser("curve", help="SF as a function of mu_r")
    sp.add_argument("--ratio", type=float, default=0.5)
    common(sp); sp.set_defaults(func=cmd_curve)

    sp = sub.add_parser("compare", help="sphere versus cylinder for several a/b")
    sp.add_argument("--mu", type=float, default=1000.0)
    sp.add_argument("--ratios", type=float, nargs="+", default=[0.8, 0.3])
    common(sp); sp.set_defaults(func=cmd_compare)

    sp = sub.add_parser("crossover", help="a/b at which both shells shield equally")
    sp.add_argument("--mus", type=float, nargs="+", default=[10, 1000, 1e4])
    sp.set_defaults(func=cmd_crossover)

    sp = sub.add_parser("table", help="numerical table of shielding factors")
    sp.add_argument("--ratio", type=float, default=0.5)
    sp.add_argument("--mus", type=float, nargs="+", default=[1, 2, 5, 10, 100, 1000])
    sp.set_defaults(func=cmd_table)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
