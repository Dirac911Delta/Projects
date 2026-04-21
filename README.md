# Projects

MATLAB (R2023) modular framework for transient thermal simulation of end-pumped Nd-doped laser glasses.

## Files
- `main.m` - user configuration and simulation driver
- `materials.m` - material database and lifetime interpolation hooks
- `pump_model.m` - super-Gaussian pump + Beer–Lambert absorption + heat source
- `heat_solver.m` - axisymmetric Crank–Nicolson implicit heat solver
- `post_processing.m` - derived metrics and plots

## Important data transparency notes
- Missing table values are stored as `NaN` (no fabrication).
- 808 nm absorption is **not available** in the provided optical table image and is handled via:
  - `settings.absorption_808_user` (preferred), or
  - concentration-proportional estimate (explicitly flagged warning, requires user-supplied reference coefficient).
- Lifetime-vs-doping data is implemented as a lookup/interpolation interface and can be user-overridden.

## Run
Open MATLAB in this folder and run:

```matlab
main
```

## Extending materials
Edit `materials.m` and populate:
- `emission_cross_section`
- `effective_bandwidth`
- `fluorescence_wavelength`
- `Nd_concentration`
- `lifetime_vs_doping.doping_wt_pct`
- `lifetime_vs_doping.lifetime_s`
- `absorption_808` (if measured)
- `thermal_conductivity`, `specific_heat`, `density`

Keep units in SI where applicable.
