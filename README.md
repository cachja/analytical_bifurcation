# Supporting Wolfram Mathematica codes for "Minimal Explicit Spectral Analysis of Oscillatory Onset in Flow Past a Cylinder and a Sphere" paper.

These Mathematica scripts are intended as deterministic reproducibility checks for the reduced models in the paper.

Files:
  1) cylinder_stokes_2modes.wl
     Frozen annulus Stokes base flow + 2 angular modes.

  2) cylinder_corrected_3modes.wl
     First-order corrected annulus Stokes base flow + 3 angular modes.

  3) sphere_stokes_2modes.wl
     Frozen sphere Stokes base flow + 2 real m=1 poloidal modes.

  4) sphere_corrected_3modes.wl
     First-order corrected sphere base flow (computed automatically in a reduced axisymmetric Stokes space) + 3 real m=1 poloidal perturbation modes.

General notes:
  - All scripts start from explicit base-flow formulas and trial modes.
  - Reduced matrices are computed by direct numerical quadrature.
  - The corrected-base scripts compute the first correction automatically.
  - CSV files for the tracked spectral branches are exported for plotting.

Suggested workflow:
  - Run the frozen-base scripts first.
  - Then run the corrected-base scripts.
  - Compare the printed matrices and the exported branch data with the numbers quoted in the manuscript.

