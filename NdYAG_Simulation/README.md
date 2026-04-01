# Nd:YAG End-Pumped Laser Simulation Suite

**Based on:** Harrison, Forbes & Naidoo, *Optics Express* **31**(15), 24516 (2023)  
"Improving performance prediction of diode end-pumped solid-state Nd:YAG rod amplifiers by incorporating pump mode evolution"

---

## Overview

This MATLAB (2023a) simulation suite implements the rigorous 3-D analytical / iterative-Fourier beam-propagation model from Harrison et al. 2023, then extends it to two additional configurations requested by the user.

### Three scenarios

| Scenario | Configuration | Crystal |
|---|---|---|
| 1 | Single-pump contra-prop amplifier (**paper baseline**) | L=25 mm, R=2 mm, 0.5% Nd:YAG |
| 2 | Dual end-pump amplifier | L=70 mm, R=3 mm, 0.5% Nd:YAG |
| 3 | CW flat/flat laser oscillator | L=70 mm, R=3 mm, 0.5% Nd:YAG |

---

## Physics Model

### 1. Pump beam shaping — phase-only DOE (Equations 1–4, Harrison 2023)

The multimode fibre-coupled (MMFC) diode pump is modelled as a complex field using a **phase-only Gaussian-to-flat-top (FT) diffractive optical element**:

```
Φ_FT(r) = (π·ωᵢ·ω_f)/(√2·f) · ∫₀^ρ √(1 − exp(−ρ'²)) dρ'    [Eq. 1]
T_lens(r) = exp(−i·k·r²/2f)                                    [Eq. 2]
DOE₁(r)   = T_lens(r) · exp(−i·Φ_FT(r))                       [Eq. 3]
Ψ_G(r)    = √(2Pₚ/π·ωᵢ²) · exp(−r²/ωᵢ²) · DOE₁(r)           [Eq. 4]
```

The modulated Gaussian `Ψ_G` evolves into a flat-top profile of radius `ω_f` over the Fourier-lens propagation distance `f`. Co-propagation uses the conjugate DOE phase (DOE₂).

### 2. Beam propagation — Angular Spectrum Method (BPM)

The pump and seed fields are propagated through the crystal using the **split-step angular-spectrum method** (Goodman, *Introduction to Fourier Optics*):

1. Half-step thin-lens mask: absorption + gain + phase from Δn
2. Free-space propagation: `H(kₓ,kᵧ) = exp(i·k_z·dz)`
3. Second half-step thin-lens mask

This allows the dynamic evolution of the pump beam profile to be captured at each z-slice.

### 3. Spectral absorption cross-section

The effective pump absorption coefficient is computed as the **spectral overlap integral** of the diode emission spectrum with the Nd:YAG absorption lineshape:

```
σ_eff = ∫ σ_abs(λ) · S(λ) dλ / ∫ S(λ) dλ
```

A temperature-induced redshift of the diode wavelength (0.27 nm/K) is included.

### 4. Population inversion — 4-level rate equations

Steady-state solution for a 4-level system (Nd³⁺ in YAG):

```
N₂ = N_tot · Rₚ / (Rₚ + 1/τ_em + σ_em·φ_s)

where:
  Rₚ = α_eff · I_p / (N_tot · hν_p)   [pump rate per ion, s⁻¹]
  φ_s = I_s / hν_s                     [signal photon flux, m⁻²s⁻¹]
  I_sat = hν_s / (σ_em · τ_em)        [saturation intensity, W/m²]
  g = σ_em · N₂                        [power gain coefficient, m⁻¹]
```

### 5. Thermal model — Kirchhoff transform + FDM

The steady-state heat equation with **temperature-dependent thermal conductivity** `k(T) = k₀·(T/T₀)^ξ` (Aggarwal et al., Table 1 parameters: k₀=10.5 W/m/K, T₀=300 K, ξ=−0.77) is solved via the Kirchhoff transform:

```
u(T) = k₀·T₀/(1+ξ) · [(T/T₀)^(1+ξ) − (T_sink/T₀)^(1+ξ)]
∇²u = −Q(r,z)     (now linear, solved by 2-D FDM per z-slice)
```

Boundary condition: Newton's cooling at crystal surface with `h_c = 1.6×10⁴ W/m²/K` (Table 1).

### 6. Refractive index change — thermal + elasto-optic ([111]-cut)

```
Δn(r,z) = (dn/dT)·ΔT  +  Δn_SO(r,z)
```

For [111]-cut Nd:YAG (Salinas-Alvarado et al., ref [17]):
- Thermoelastic stresses σ_r, σ_φ computed from the Lamé solution for a heated cylinder
- Effective elasto-optic coefficient: `B_eff = (p₁₁ + 2p₁₂ + 4p₄₄)/3`  
  with p₁₁=−0.029, p₁₂=0.0091, p₄₄=−0.0615 (Table 1)
- Thermal lens focal length from parabolic fit to accumulated phase `φ_TL = k₀∫Δn dz`

### 7. CW cavity (Scenario 3) — Fox-Li iteration

Round-trip BPM propagation iterated until convergence:
- Forward pass (A→B) with gain + thermal phase → apply OC (√R_OC)
- Backward pass (B→A) with gain + thermal phase → apply HR (√R_HR)
- Threshold: `exp(2∫g dz) · R_OC · R_HR · (1−losses) = 1`

---

## Crystal Parameters (Table 1, Harrison 2023)

| Parameter | Symbol | Value | Unit |
|---|---|---|---|
| Crystal length (paper) | L | 25 | mm |
| Crystal radius (paper) | R | 2 | mm |
| Atomic dopant | – | 0.5 | at.% |
| Active ion density | N_tot | 0.69×10²⁰ | cm⁻³ |
| Refractive index | n₀ | 1.82 | – |
| Fluorescence lifetime | τ_em | 240 | µs |
| Emission cross-section | σ_em | 2.8×10⁻¹⁹ | cm² |
| Absorption XS (peak) | σ_abs^peak | 5×10⁻²⁰ | cm² |
| Absorption XS (dynamic) | σ_abs^p | 5.9–9.1×10⁻²¹ | cm² |
| Thermal conductivity | k₀ | 10.5 | W/m/K |
| Heat-sink temperature | T_sink | 294 | K |
| Thermal conductance | h_c | 1.6×10⁴ | W/m²/K |
| Young's modulus | M | 3×10¹⁰ | Pa |
| Poisson's ratio | ν | 0.25 | – |
| Elasto-optic p₁₁ | p₁₁ | −0.029 | – |
| Elasto-optic p₁₂ | p₁₂ | 0.0091 | – |
| Elasto-optic p₄₄ | p₄₄ | −0.0615 | – |

---

## File Structure

```
NdYAG_Simulation/
├── run_simulation.m               ← Top-level entry point (run this)
│
├── NdYAG_params.m                 ← All physical parameters (Table 1)
├── compute_DOE_phase.m            ← DOE pump beam shaping (Eq. 1–4)
├── compute_absorption_xsec.m      ← Spectral absorption cross-section
├── bpm_propagate.m                ← Angular-spectrum BPM (split-step)
├── solve_heat_2d.m                ← 3-D thermal model (Kirchhoff + FDM)
├── compute_refractive_index_change.m  ← Δn thermal + elasto-optic [111]
├── compute_population_inversion.m ← 4-level steady-state rate equations
│
├── single_pump_amplifier.m        ← Scenario 1: paper baseline
├── dual_pump_amplifier.m          ← Scenario 2: dual end-pump
├── CW_flat_flat_laser.m           ← Scenario 3: CW flat/flat cavity
└── README.md                      ← This file
```

---

## Quick Start

```matlab
% Run all three scenarios
run_simulation

% Or run individually:
res1 = single_pump_amplifier();           % paper baseline
res2 = dual_pump_amplifier();             % dual pump, user crystal
res3 = CW_flat_flat_laser();              % CW flat/flat cavity

% Use the user's larger crystal for Scenario 1:
res1u = single_pump_amplifier([], [], 'contra', 'user');

% Sweep a custom pump power range:
res1c = single_pump_amplifier(0:5:38, 0.1, 'contra', 'paper');
```

---

## Requirements

- MATLAB R2023a or later  
- No additional toolboxes required (all built-in functions only)

---

## References

1. J. Harrison, A. Forbes, D. Naidoo, *Opt. Express* **31**(15), 24516 (2023)  
2. W. Koechner, *Solid-State Laser Engineering*, 6th ed. (Springer, 2006)  
3. D. Salinas-Alvarado et al., *J. Appl. Phys.* **111**, 013112 (2012)  
4. R. L. Aggarwal et al., *J. Appl. Phys.* **98**, 103514 (2005)  
5. J. W. Goodman, *Introduction to Fourier Optics*, 3rd ed. (Roberts, 2005)
