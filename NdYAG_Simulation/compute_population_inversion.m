function [N2, gain_2d] = compute_population_inversion(p, I_pump_2d, I_signal_2d, alpha_eff)
% COMPUTE_POPULATION_INVERSION  Steady-state upper-level population N₂ and
%   the resulting gain coefficient g for a 4-level Nd:YAG system.
%
% 4-Level rate equations (steady state, dN₂/dt = 0):
% -------------------------------------------------------
%   dN₂/dt = Rp·(N_tot − N₂) − N₂/τ_em − σ_em·φ_s·N₂ = 0
%
% where:
%   Rp    = σ_abs · I_p / hν_p  = α_eff · I_p / (N_tot · hν_p)  [s⁻¹]
%           pump absorption rate per active ion
%   φ_s   = I_s / hν_s                                           [m⁻² s⁻¹]
%           signal photon flux
%   τ_em  = fluorescence lifetime (240 µs, Table 1)
%   I_sat = hν_s / (σ_em · τ_em)                                 [W/m²]
%           signal saturation intensity
%
% Solving for N₂:
%   N₂ = N_tot · Rp · τ_em / (Rp·τ_em + 1 + σ_em·φ_s·τ_em)
%      = N_tot · Rp / (Rp + 1/τ_em + I_s/I_sat · 1/τ_em)
%
% Gain coefficient (4-level, no ground-state absorption at 1064 nm):
%   g(r,z) = σ_em · N₂(r,z)   [m⁻¹]
%
% Reference: Section 3.2 rate equations, Harrison et al. 2023.
%
% INPUTS
%   p           - parameter struct (NdYAG_params)
%   I_pump_2d   - [N×N]  pump intensity at this z-slice        [W/m²]
%   I_signal_2d - [N×N]  signal intensity at this z-slice      [W/m²]
%                        (use zeros for small-signal / no signal)
%   alpha_eff   - effective pump absorption coefficient         [m⁻¹]
%
% OUTPUTS
%   N2      - [N×N] upper-level population density              [m⁻³]
%   gain_2d - [N×N] power gain coefficient g = σ_em · N₂       [m⁻¹]

N_tot    = p.N_tot;
sigma_em = p.sigma_em;
tau_em   = p.tau_em;
hnu_p    = p.hnu_p;
hnu_s    = p.hnu_s;
I_sat    = p.I_sat;

% ---- Pump absorption rate per ion  Rp = α·I_p / (N_tot · hν_p) [s⁻¹] -
Rp = (alpha_eff .* I_pump_2d) ./ (N_tot .* hnu_p);

% ---- Signal saturation term  Is/I_sat / τ_em  (= σ_em · φ_s) -----------
sat_term = I_signal_2d ./ (I_sat .* tau_em);   % = σ_em · φ_s  [s⁻¹]

% ---- Steady-state N₂ ---------------------------------------------------
%   N₂ = N_tot · Rp · τ_em / (Rp·τ_em + 1 + σ_em·φ_s·τ_em)
%   Divide numerator and denominator by τ_em:
%   N₂ = N_tot · Rp / (Rp + 1/τ_em + σ_em·φ_s)

denom = Rp + (1/tau_em) + sat_term;          % [s⁻¹]
N2    = N_tot .* Rp ./ denom;               % [m⁻³]
N2    = max(N2, 0);                          % physical floor

% ---- Small-signal power gain coefficient  g = σ_em · N₂  [m⁻¹] --------
gain_2d = sigma_em .* N2;

end
