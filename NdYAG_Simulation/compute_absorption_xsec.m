function alpha_eff = compute_absorption_xsec(p, diode_T_offset)
% COMPUTE_ABSORPTION_XSEC  Compute the spectrally-weighted effective pump
%   absorption coefficient α_eff [m⁻¹] for the MMFC diode pump source.
%
% The MMFC diode emission spectrum is approximated as a Gaussian centred on
% λ_p with FWHM = 2·Δλ (from Table 1).  The Nd:YAG absorption lineshape at
% 808 nm is also Gaussian with a narrower FWHM.  The effective absorption
% cross-section is the spectral overlap integral:
%
%   σ_eff = ∫ σ_abs(λ) · S(λ) dλ  /  ∫ S(λ) dλ
%
% where S(λ) is the normalised diode emission power spectral density.
%
% A temperature-induced redshift of the diode wavelength is included:
%   δλ/δT ≈ +0.27 nm/K  (typical for AlGaAs diodes at 808 nm)
%
% Reference: Section 3.2, Harrison et al., Opt. Express 31(15), 24516 (2023)
%            Table 1 for spectral parameters.
%
% INPUTS
%   p               - parameter struct (NdYAG_params)
%   diode_T_offset  - diode junction temperature rise above nominal [K]
%                     (default 0 → no shift)
%
% OUTPUT
%   alpha_eff       - effective absorption coefficient [m⁻¹]
%                     α_eff = σ_eff · N_tot

if nargin < 2 || isempty(diode_T_offset), diode_T_offset = 0; end

% ---- Spectral grid (fine λ-grid around pump band) -----------------------
lambda_c  = p.lambda_p + 0.27e-9 * diode_T_offset; % shifted centre [m]
dlambda   = p.delta_lambda_p;                        % diode half-width [m]
Nlam      = 4096;
lam_min   = lambda_c - 5 * dlambda;
lam_max   = lambda_c + 5 * dlambda;
lam_vec   = linspace(lam_min, lam_max, Nlam);        % [m]

% ---- Diode emission spectrum S(λ) – Gaussian (Section 2 / Table 1) -----
% FWHM = 2·δλ = 2·1.85 nm → σ_diode = FWHM/(2√(2ln2))
sigma_diode = (2 * dlambda) / (2 * sqrt(2 * log(2)));
S = exp(-0.5 * ((lam_vec - lambda_c) / sigma_diode).^2);  % normalised

% ---- Nd:YAG absorption lineshape σ_abs(λ) at 808 nm -------------------
% The absorption peak at 808 nm has a FWHM of ≈ 1.5 nm (homogeneous +
% inhomogeneous broadening in YAG at room temperature).
lam_peak_abs  = 808e-9;                   % absorption line centre [m]
FWHM_abs      = 1.5e-9;                   % absorption FWHM [m]
sigma_abs_lw  = FWHM_abs / (2 * sqrt(2 * log(2)));
% Lorentzian × Gaussian (Voigt) approximated as Gaussian for simplicity
sigma_abs_spectrum = p.sigma_abs_peak .* ...
    exp(-0.5 * ((lam_vec - lam_peak_abs) / sigma_abs_lw).^2);

% ---- Spectrally-weighted effective cross-section -----------------------
sigma_eff = trapz(lam_vec, sigma_abs_spectrum .* S) / trapz(lam_vec, S);

% Clamp to the measured dynamic range from Table 1
sigma_eff = max(p.sigma_abs_dyn_lo, min(p.sigma_abs_dyn_hi, sigma_eff));

% ---- Effective absorption coefficient ----------------------------------
alpha_eff = sigma_eff * p.N_tot;    % [m⁻¹]

end
