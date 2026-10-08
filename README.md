# Magnetic shielding by permeable shells: exact fields, interactive notebook and classroom activities

Companion code to the article *Magnetic shielding as a motivational tool in
teaching classical electromagnetic theory* (submitted to the European Journal
of Physics).

The code computes and plots the exact static field of a permeable spherical
shell and of an infinitely long cylindrical shell in a uniform applied field,
treats both geometries in a single calculation, and reproduces every figure of
the article and of its Supplementary Material II. All scripts run offline and
depend only on `numpy` and `matplotlib` (plus `ipywidgets` for the notebook).

## The physics in one paragraph

A linear, isotropic and homogeneous shell of relative permeability `mu_r`,
inner radius `a` and outer radius `b` is placed in vacuum in a static uniform
field `H0` (so that `B0 = mu0 H0` far away). The field inside the cavity is
exactly uniform and parallel to `H0`, and the shielding factor
`SF = H0/H_int` of both shells takes one form,

    SF = 1 + [(d - 1)/d^2] * [1 - (a/b)^d] * (mu_r - 1)^2 / mu_r,

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
| `verify_revision_claims.py` | Checks every quantitative statement of the article and of its supplements, and the consistency of the code (86 checks). |
| `verify_cavity_uniformity.py` | Independent check, by a finite-volume solution, that the cavity field is exactly uniform. |
| `build_notebook.py` | Maintainer tool: rebuilds the notebook from `shielding.py`. |
| `requirements.txt` | Python packages. |

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
condition, every number quoted in the article, the continuity of the
tangential component of `H` and of the normal component of `B` for the plotted
fields, and the agreement between `shielding.py`, an independent
implementation and the notebook. The second solves the radial problem by finite
volumes and shows that the cavity field is uniform to rounding error.

## Conventions and idealizations

The sphere is shown in a meridian plane and the cylinder in its transverse
plane, with the applied field pointing to the right. Field strengths are given
in units of `H0` and flux densities in units of `mu0 H0`; the white curves of
the maps are streamlines of `B`, whose spacing carries no information about the
field strength. The model assumes a linear, isotropic and homogeneous material
with a field-independent permeability, a static uniform applied field, an infinitely long
cylinder and closed shells without holes or seams; saturation, hysteresis,
finite length and apertures are outside its scope.

## History

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

Parts of the code were refactored with the assistance of Claude (Anthropic).
All code was reviewed and tested by the authors, who take full responsibility
for it.

## Citation

If you use this code, please cite the article (reference to be added upon
publication) and the archived version of this repository (DOI to be added).
