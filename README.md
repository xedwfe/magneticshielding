# Magnetic shielding by permeable shells: exact fields, interactive notebook and classroom activities

Companion code to the article *Magnetic shielding as a context for teaching
magnetostatics in matter: spherical and cylindrical shells with open code*
(submitted to the European Journal of Physics).

The code computes and plots the exact static field of a permeable spherical
shell and of an infinitely long cylindrical shell in a uniform applied field,
treats both geometries in a single calculation, and reproduces every figure of
the article and the code-generated figures of its two supplements. All scripts
run offline and depend only on `numpy` and `matplotlib` (plus `ipywidgets` for
the notebook). The code is released under the MIT License (see `LICENSE`).

## The physics in one paragraph

A linear, isotropic and homogeneous shell of relative permeability `mu_r`,
inner radius `a` and outer radius `b` is placed in vacuum in a static uniform
field `H0` (so that `B0 = mu0 H0` far away). The field inside the cavity is
exactly uniform and parallel to `H0`, and the shielding factor
`SF = H0/H_int` of both shells takes one form,

    SF = 1 + C_d * (1 - k_d) * (mu_r - 1)^2 / mu_r,   C_d = (d - 1)/d^2,   k_d = (a/b)^d,

with `d = 3` for the sphere and `d = 2` for the cylinder in a transverse field.
The sphere is the better shield only for thin shells: the two shielding factors
are equal at `a/b = (1 + sqrt(33))/16 = 0.4215...` for every value of `mu_r`,
and the cylinder is better for thicker shells. Supplementary Material I of the
article derives these results.

## Quick start

**Without installation.** Open the notebook directly in a web browser and run
all cells; it is self-contained.

* Binder (no account needed):
  https://mybinder.org/v2/gh/xedwfe/magneticshielding/HEAD?labpath=magnetic_shielding_interactive.ipynb
* Google Colab (requires a Google account):
  https://colab.research.google.com/github/xedwfe/magneticshielding/blob/main/magnetic_shielding_interactive.ipynb

**Locally.**

    pip install -r requirements.txt jupyter
    jupyter notebook magnetic_shielding_interactive.ipynb

## Files

| File | Purpose |
|---|---|
| `magnetic_shielding_interactive.ipynb` | Interactive notebook: field maps of both shells with sliders for `mu_r`, `a/b`, `H0`, `b`, the window half-width `L` and the grid resolution `N`; sections on the shielding-factor curve, on the comparison between sphere and cylinder, and on `B` versus `H`. |
| `shielding.py` | The exact solution and derived quantities (single source of the physics). |
| `paper_figures.py` | Regenerates every figure of the article and the code-generated figures of Supplementary Materials I and II as vector PDF. |
| `explore_shielding.py` | Command-line exploration without Jupyter (field maps, tables, curves, crossing point). |
| `verify_revision_claims.py` | Checks every quantitative statement of the article and of its supplements, and the consistency of the code (96 checks). |
| `verify_cavity_uniformity.py` | Independent check, by a finite-volume solution, that the cavity field is exactly uniform. |
| `build_notebook.py` | Maintainer tool: rebuilds the notebook from `shielding.py`. |
| `requirements.txt` | Python packages. |
| `LICENSE` | MIT License. |

## Reproducing the figures

    python paper_figures.py

| Output | Figure |
|---|---|
| `fig1_geometry.pdf` | Article, figure 1 |
| `fig2_field_maps.pdf` | Article, figure 2 |
| `fig3_shielding_factor.pdf` | Article, figure 3 |
| `sm1_cylinder_section.pdf` | Supplementary Material I, figure 4 |
| `sm1_sphere_section.pdf` | Supplementary Material I, figure 5 |
| `sm2_fig2_mu_sweep.pdf` | Supplementary Material II, figure 2 |
| `sm2_fig3_geometry.pdf` | Supplementary Material II, figure 3 |
| `sm2_fig4_B_vs_H.pdf` | Supplementary Material II, figure 4 |

Figure 1 of Supplementary Material II is a screenshot of the notebook, and
figures 1 to 3 of Supplementary Material I are hand-drawn illustrations.

The figures of the article were produced with Python 3.13, `numpy` 2.5 and
`matplotlib` 3.11. Other versions give the same fields and the same numbers,
but the placement of the streamlines, which carries no physical information,
may differ slightly between versions of `matplotlib`.

## Exploring without Jupyter

    python explore_shielding.py maps --mu 1000 --ratio 0.5                 # field maps of both shells
    python explore_shielding.py maps --mu 1000 --ratio 0.5 --quantity H    # |H| instead of |B|
    python explore_shielding.py table --ratio 0.5 --mus 1 2 10 100 1000    # numerical table
    python explore_shielding.py curve --ratio 0.5                          # SF as a function of mu_r
    python explore_shielding.py compare --mu 1000 --ratios 0.8 0.3         # sphere versus cylinder
    python explore_shielding.py crossover --mus 10 1000 10000              # crossing point

Each command prints its numbers; the plotting commands also save a PNG file
(`--out` chooses the name, `--show` opens a window).

## Verification

    python verify_revision_claims.py
    python verify_cavity_uniformity.py

The first script checks the closed forms against direct solutions of the
interface conditions, the crossing point, the limits and their reading in
terms of demagnetizing factors, the shielding of higher harmonics, the cloaking
condition, every number quoted in the article and in its supplements, the continuity of the
tangential component of `H` and of the normal component of `B` for the plotted
fields, and the agreement between `shielding.py`, an independent
implementation and the notebook. The second lists the determinants of the
interface conditions for the harmonics absent from a uniform applied field
(`l = 0, 2, 3, ...` for the sphere and `nu = 0, 2, 3, ...` for the cylinder),
solves the radial problem by finite volumes and shows that the cavity field is
uniform to rounding error.

## Notation

The notation is that of the article and its supplements. The shell has inner
radius `a`, outer radius `b`, thickness `t = b - a` and relative permeability
`mu_r`; its wall is the region `a < r < b`, where `r` is the distance from the
centre of the sphere or from the axis of the cylinder (`rho` in the article).
`H_int` is the strength of the uniform cavity field and `SF = H0/H_int` the
shielding factor; `C_d = (d - 1)/d^2`, `k_d = (a/b)^d`, and `N = 1/d` is the
demagnetizing factor of the corresponding solid body. The angle `theta` is
measured from the direction of `H0`; the poles of each surface are its points
with `theta = 0` and `theta = pi`, and its equator is formed by the points with
`theta = pi/2`. The potentials in the cavity, in the shell and outside are
`phi_in`, `phi_shell` and `phi_out`.

## Conventions and idealizations

The sphere is shown in a meridional plane and the cylinder in its transverse
plane, with the applied field pointing to the right. In the figures of the
article (`paper_figures.py`, `explore_shielding.py`) this axis is labelled `z`
for the sphere, its symmetry axis, and `x` for the cylinder; the notebook
labels it `x` in both panels. Field strengths are given in units of `H0` and
flux densities in units of `mu0 H0`; the white curves of the maps are
streamlines of `B`, whose spacing carries no information about the field
strength. The model assumes a linear, isotropic and homogeneous material
with a field-independent permeability, a static uniform applied field, an infinitely long
cylinder and closed shells without holes or seams; saturation, hysteresis,
finite length and apertures are outside its scope.

## History

**Version 3.3 (2026, submitted revision).**
`LICENSE` (MIT) and `.gitignore` added. The boxes of the notebook maps write
`H_int/H0` in scientific notation (for example 5.13 x 10^-3), as in the figure
of Supplementary Material II that shows the notebook. This file states the
software versions used for the figures of the article.

**Version 3.2 (2026, notation aligned with the final text).**
Article title updated. `shielding.py`, `paper_figures.py`, the notebook and this
file use the notation of the article and its supplements: meridional plane,
potentials `phi_in`, `phi_shell` and `phi_out`, the angle `theta` measured from
`H0`, poles and equator, the wall `a < r < b`, and a superconducting layer of
outer radius `a` in the magnetic cloak; `demagnetising_factor` is renamed
`demagnetizing_factor` (the old name is kept as an alias). The dotted line of
article figure 3(a) is labelled `SF = C_3(1 - k_3) mu_r`, as in its caption, and
the box of the notebook maps writes `H_int` in roman type. `verify_revision_claims.py`
covers the statements added in this revision (96 checks), including the
harmonics `l = 0` of the sphere and `nu = 0` of the cylinder, which
`verify_cavity_uniformity.py` now also lists.

**Version 3.1 (2026, final revision of the article).**
`verify_revision_claims.py` extended to every quantitative statement of the
article and of both supplements (86 checks). `paper_figures.py` adds the polar
coordinates to the two sections of Supplementary Material I and now also
produces its figure 4. Direct Binder and Colab links added to this file.

**Version 3.0 (2026, revision of the article).**
`shielding.py` introduced, with one derivation for both geometries.
`paper_figures.py` produces all figures as vector PDF at printed size, with
`|B|` in colour and streamlines of `B`, and adds the ratio of the two shielding
factors as a function of `a/b`. `explore_shielding.py` replaces
`magnetic_shieldingV3_corrected.py` and `magnetic_shielding_S_graf_corrected.py`.
The notebook was rewritten, keeping the same controls, with new sections on
the shielding-factor curve, on the geometry comparison and on `B` versus `H`.
Corrections: the sphere is not always the better shield (the ordering reverses
at `a/b = 0.42`); the overestimate of the approximation `SF = A mu_r + 1` is 19%
at `mu_r = 1`, at most 27% (near `mu_r = 2`) and 14% at `mu_r = 10`; the ranges of
permeability of soft ferrites, silicon steels and nickel-iron alloys.

**Version 2.0.** Exact shielding factors instead of the approximation `+1`;
correction of the exterior coefficient of the spherical shell; interior field
lines shown; corrected material annotations and reference to Hoburg (1995).

## Use of AI tools

During the preparation of this work the authors used Claude (Anthropic) to
assist in refactoring the simulation code, in drafting and editing the text of
the article and of its supplements, and in cross-checking algebraic results
numerically. All derivations, code and results were reviewed, tested and
verified by the authors, who take full responsibility for them.

## License

MIT License; see `LICENSE`.

## Citation

If you use this code, please cite the article (reference to be added upon
publication) and the archived version of this repository (DOI to be added).
