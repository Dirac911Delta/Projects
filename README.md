# Nd Glass Thermal Simulation – MATLAB R2023

Modular transient thermal simulation for end-pumped Nd-doped laser glasses.

## Files

| File | Purpose |
|---|---|
| `main.m` | User controls, simulation driver, console summary |
| `materials.m` | Complete material database (all three manufacturer tables) |
| `pump_model.m` | Super-Gaussian pump, Beer–Lambert absorption, quantum-defect heating |
| `heat_solver.m` | Axisymmetric Crank–Nicolson implicit thermal solver |
| `post_processing.m` | Thermal lens, stress, safety checks, gain, plots |

## Quick start

```matlab
main      % runs with default settings (N31, 50 W, 1.0 wt% doping)
```

## How material data is used

Every field extracted from the three manufacturer tables is wired into physics:

| Field | Used in |
|---|---|
| `thermal_conductivity` | heat_solver – thermal diffusivity |
| `specific_heat` | heat_solver – thermal diffusivity |
| `density` | heat_solver – thermal diffusivity |
| `fluorescence_wavelength` | pump_model – quantum-defect heating fraction η = 1−λ_p/λ_fl |
| `refractive_index_1053nm` | pump_model – optional Fresnel face-reflection loss |
| `emission_cross_section` | pump_model – saturation intensity; post_processing – gain |
| `dn_dT` | post_processing – thermal lens (fallback when W₀ not available) |
| `thermal_opl_coeff` (W₀) | post_processing – full thermal lens focal power via parabolic fit |
| `thermal_expansion_coeff` | post_processing – Timoshenko thermal stress; end-face bulging |
| `youngs_modulus_Pa` | post_processing – radial and hoop thermal stress |
| `poissons_ratio` | post_processing – thermal stress |
| `fracture_toughness_MPa_sqrtm` | post_processing – fracture risk reference |
| `transition_temperature_C` | post_processing – peak-temperature safety warning |
| `softening_temperature_C` | post_processing – peak-temperature safety warning |
| `Nd_concentration` | post_processing – small-signal gain coefficient g₀ = σ·N |
| `absorption_coeff_1053nm_max` | post_processing – internal transmission at 1053 nm |
| `nonlinear_n2_esu` | pump_model – approximate B-integral |
| `lifetime_vs_doping` | main / pump_model – saturation intensity, lifetime summary |
| `effective_bandwidth` | console summary (gain bandwidth) |
| `abbe_value` | stored; relevant for chromatic dispersion calculations |
| `wt_pct_nominal` | stored; reference nominal doping from table |

## Critical missing data note

**808 nm absorption** is not present in any manufacturer table.  
Set `settings.absorption_808_user` to a measured value, or enable
`settings.estimate_absorption_from_doping = true` and provide a reference
coefficient and concentration. Any estimation is explicitly flagged with a warning.

## Extending to a new material

Add an entry to `build_database()` in `materials.m` following the same
index-parallel vector layout. Set any unknown field to `NaN`. Populate
`lifetime_tables` with `{[doping_wt_pct], [lifetime_us]}` pairs.

## Physical assumptions flagged in code

- Quantum-defect fraction `η_heat = 1 − λ_pump/λ_fl` (Stokes limit; excludes non-radiative QE)
- Fresnel loss uses `n(1053 nm)` as proxy for `n(808 nm)` (≈ 0.5–1% error for phosphate glass)
- B-integral conversion `n₂[m²/W] ≈ n₂[1e−13 e.s.u.] × 1e−13 × 4.19e−7/n₀`  
  Reference: Hellwarth R.W. (1977) *Progress in Quantum Electronics* **5**, 1–68.
- Thermal lens parabolic fit over inner 1/3 of aperture; see Innocenzi et al., *J. Appl. Phys.* **75**, 4991 (1994)
- Thermal stress uses plane-stress Timoshenko cylinder formula; valid for free-ended rods
- Generic thermal fallbacks (k = 0.80 W/(m·K), cp = 800 J/(kg·K), ρ = 2700 kg/m³) used only when table value is NaN
