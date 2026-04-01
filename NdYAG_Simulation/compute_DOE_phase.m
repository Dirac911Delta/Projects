function [E_pump, alpha_eff] = compute_DOE_phase(p, P_pump, propagation_dir, diode_T)
% COMPUTE_DOE_PHASE  Generate the pump field at the crystal input face using
%   the phase-only Gaussian-to-flat-top (FT) diffractive optical element
%   described in Harrison et al., Opt. Express 31(15), 24516 (2023).
%
% The MMFC diode pump beam is modelled as a complex field that starts as a
% Gaussian at the fibre exit and is phase-modulated by a phase-only DOE so
% that it transforms into a flat-top (FT) profile of radius ω_f at the
% far surface of the crystal during free-space (Fourier-lens) propagation.
%
% === Equations from the paper ===========================================
%
% Phase-only Gaussian-to-FT transmission function  [Eq. 1, Harrison 2023]:
%
%   Φ_FT(r) = (π·ω_i·ω_f) / (√2·f) · ∫₀^ρ √(1 − exp(−ρ'²)) dρ'
%
%   where  ρ  = √2·r / ω_i ,
%          ω_i = initial Gaussian beam radius (at DOE plane),
%          ω_f = desired FT beam radius (at back focal plane of lens f),
%          f   = focal length of Fourier lens.
%
% Fourier lens transmittance  [Eq. 2]:
%   T_lens(r) = exp(−i·k·r² / (2f))
%
% DOE for contra-propagating configuration  [Eq. 3]:
%   DOE₁(r) = T_lens(r) · exp(−i·Φ_FT(r))
%
% Initial pump field (contra-prop) [Eq. 4]:
%   Ψ_G(r) = √(2·P_p / (π·ω_i²)) · exp(−r²/ω_i²) · DOE₁(r)
%
% For co-propagating configuration, the pump enters as a FT beam at the
% same face as the seed.  The DOE₂ is the complex conjugate phase applied
% on the FT profile so that it refocuses to ω_i at z = L  [Section 3.1.2]:
%   DOE₂(r) = T_lens*(r) · exp(+i·Φ_FT(r))
%   Ψ_FT(r) = √(P_p / (π·ω_f²)) · circ(r/ω_f) · DOE₂(r)
%
% ========================================================================
%
% INPUTS
%   p              - parameter struct from NdYAG_params
%   P_pump         - pump power                             [W]
%   propagation_dir- 'contra' (pump enters at face A, seed at B) [default]
%                    'co'     (pump and seed both enter at A)
%   diode_T        - diode temperature offset for wavelength shift [K]
%                    (default 0 → use p.lambda_p)
%
% OUTPUTS
%   E_pump         - [N×N] complex pump field at crystal input face
%                    |E|² = intensity [W/m²]
%   alpha_eff      - effective (spectrally-weighted) absorption coeff [m⁻¹]
%                    computed via compute_absorption_xsec
%
% --------------------------------------------------------------------------

if nargin < 3 || isempty(propagation_dir), propagation_dir = 'contra'; end
if nargin < 4 || isempty(diode_T),         diode_T = 0;                end

N    = p.N_grid;
R    = p.R_grid;          % radial coordinate map  [m]
k    = p.k_p;             % pump wavenumber in crystal
omega_i = p.omega_i;      % Gaussian beam radius at DOE plane  [m]
omega_f = p.omega_f;      % desired FT beam radius             [m]
f       = p.f_lens_DOE;   % Fourier lens focal length          [m]

% ---- Effective absorption cross-section (spectral integral) -------------
alpha_eff = compute_absorption_xsec(p, diode_T);

% ---- Φ_FT(r)  [Eq. 1] ---------------------------------------------------
% Normalised radial variable ρ = √2·r/ω_i
rho = sqrt(2) .* R ./ omega_i;        % [N×N]

% Numerical integration  ∫₀^ρ √(1 - exp(-ρ'²)) dρ'
% using cumulative trapezoidal rule on the 1-D radial grid for efficiency,
% then map back to 2-D via interpolation.
r_1d   = linspace(0, max(R(:)), 4096);
rho_1d = sqrt(2) .* r_1d ./ omega_i;
integrand_1d = sqrt(max(1 - exp(-rho_1d.^2), 0));   % ≥0 for safety
I_1d   = cumtrapz(rho_1d, integrand_1d);             % cumulative integral over ρ'

Phi_FT_2d = (pi .* omega_i .* omega_f) ./ (sqrt(2) .* f) .* ...
            interp1(r_1d, I_1d, R, 'linear', 0);     % interpolate to 2-D

% ---- Fourier lens transmittance  T_lens(r)  [Eq. 2] ---------------------
T_lens = exp(-1i .* k .* R.^2 ./ (2 .* f));         % [N×N]

% ---- Build field at crystal input face ----------------------------------
switch lower(propagation_dir)

    case 'contra'
        % Eq. 3: DOE₁ = T_lens · exp(−i·Φ_FT)
        DOE1 = T_lens .* exp(-1i .* Phi_FT_2d);

        % Eq. 4: Ψ_G = √(2P/(π·ω_i²)) · exp(−r²/ω_i²) · DOE₁
        % Gaussian amplitude normalised so ∫|E|² dA = P_pump
        A_G  = sqrt(2 .* P_pump ./ (pi .* omega_i^2));
        E_pump = A_G .* exp(-R.^2 ./ omega_i^2) .* DOE1;

    case 'co'
        % Co-propagating: pump enters as flat-top at face A (z = 0).
        % DOE₂ is the phase-conjugate so the beam refocuses to ω_i at z = L.
        % DOE₂ = T_lens* · exp(+i·Φ_FT)  [Section 3.1.2]
        DOE2 = conj(T_lens) .* exp(+1i .* Phi_FT_2d);

        % FT amplitude (super-Gaussian order 10 ≈ flat-top within ω_f)
        A_FT = sqrt(P_pump ./ (pi .* omega_f^2));
        E_flat = A_FT .* exp(-(R ./ omega_f).^10);   % soft aperture

        % Re-normalise to exact power
        I_flat = abs(E_flat).^2;
        norm_factor = sqrt(P_pump / (sum(I_flat(:)) * p.dx^2));
        E_pump = E_flat .* norm_factor .* DOE2;

    otherwise
        error('compute_DOE_phase: propagation_dir must be "contra" or "co".');
end

% ---- Apply crystal aperture (crop field to crystal radius) --------------
E_pump = E_pump .* p.mask_crystal;

% ---- Renormalise to exact pump power after aperture --------------------
I_pump = abs(E_pump).^2;
current_power = sum(I_pump(:)) * p.dx^2;
if current_power > 0
    E_pump = E_pump .* sqrt(P_pump / current_power);
end

end
