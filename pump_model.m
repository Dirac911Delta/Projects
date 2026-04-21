function pump = pump_model(material, settings)
% pump_model.m
% End-pumped 808 nm model with physically grounded heat source.
%
% Physics implemented:
%   I(r)     = I0 * exp(-2*(r/w)^n)          super-Gaussian radial profile
%   I(r,z)   = I(r) * exp(-alpha_808 * z)     Beer-Lambert axial attenuation
%   eta_heat = 1 - lambda_pump/lambda_fl       quantum-defect heating fraction
%   Q(r,z)   = eta_heat * alpha_808 * I(r,z)  volumetric heat source [W/m^3]
%
% Optional: Fresnel reflection loss at both crystal faces uses n at 1053 nm
%           as an approximation for n at 808 nm (flagged in output).
%
% Data used from material struct:
%   fluorescence_wavelength   -> quantum-defect fraction
%   refractive_index_1053nm   -> Fresnel face-reflection loss (approx. for 808 nm)
%   emission_cross_section    -> saturation intensity
%   absorption_808            -> primary alpha source if available
%   nonlinear_n2_esu          -> B-integral estimate

PLANCK_H  = 6.62607015e-34;   % [J·s]
SPEED_C   = 2.99792458e8;     % [m/s]

nr = settings.nr;
nz = settings.nz;
R  = settings.crystal_diameter_m / 2;
L  = settings.crystal_thickness_m;
w  = settings.pump_radius_m;
P  = settings.pump_power_W;
order = settings.supergaussian_order;
lambda_pump = settings.pump_wavelength_m;      % [m]

r = linspace(0, R, nr).';
z = linspace(0, L, nz);

% -------------------------------------------------------------------------
% 1. Fresnel reflection loss at both AR-uncoated crystal faces (optional)
% -------------------------------------------------------------------------
% Uses n at 1053 nm as approximation for n at 808 nm since 808 nm index
% is not provided in the tables. Flag this assumption.
T_fresnel = 1.0;   % transmittance factor (both faces combined)
if settings.apply_fresnel_reflection_loss
    n_approx = material.refractive_index_1053nm;
    if ~isnan(n_approx)
        R_face   = ((n_approx - 1) / (n_approx + 1))^2;   % single-face reflectance
        T_fresnel = (1 - R_face)^2;                         % both faces (in, out)
        warning(['ASSUMPTION: Fresnel loss uses n(1053 nm) = %.4f as proxy for n(808 nm). ', ...
                 'Effective pump transmission factor: %.4f'], n_approx, T_fresnel);
    else
        warning('Fresnel reflection requested but refractive_index_1053nm is NaN. Skipping correction.');
    end
end
P_eff = P * T_fresnel;    % [W] effective pump power entering crystal

% -------------------------------------------------------------------------
% 2. Resolve alpha_808 [1/m]
% -------------------------------------------------------------------------
if ~isnan(settings.absorption_808_user)
    alpha = settings.absorption_808_user;
elseif ~isnan(material.absorption_808)
    alpha = material.absorption_808;
elseif settings.estimate_absorption_from_doping
    missing = {};
    if isnan(settings.absorption_808_ref),             missing{end+1} = 'absorption_808_ref'; end %#ok<AGROW>
    if isnan(settings.Nd_concentration_ref_wt_pct),    missing{end+1} = 'Nd_concentration_ref_wt_pct'; end %#ok<AGROW>
    if isnan(settings.Nd_concentration_wt_pct),        missing{end+1} = 'Nd_concentration_wt_pct'; end %#ok<AGROW>
    if ~isempty(missing)
        error('Cannot estimate alpha_808. Missing fields: %s', strjoin(missing, ', '));
    end
    % Assumption: alpha_808 scales linearly with Nd3+ concentration (wt%)
    alpha = settings.absorption_808_ref * ...
            (settings.Nd_concentration_wt_pct / settings.Nd_concentration_ref_wt_pct);
    warning(['ASSUMPTION: alpha_808 estimated by linear scaling with doping. ', ...
             'Provide settings.absorption_808_user for a measured value.']);
else
    error('alpha_808 unavailable. Set settings.absorption_808_user or enable estimate_absorption_from_doping.');
end

% -------------------------------------------------------------------------
% 3. Quantum-defect heating fraction
%    eta_heat = 1 - lambda_pump / lambda_fluorescence
%    All absorbed pump energy not emitted as a fluorescence photon (at the
%    longer wavelength) becomes heat via multi-phonon relaxation.
% -------------------------------------------------------------------------
eta_heat = 1.0;   % default: all absorbed energy becomes heat
if settings.apply_quantum_defect
    lambda_fl = material.fluorescence_wavelength;   % [m]
    if ~isnan(lambda_fl) && lambda_fl > lambda_pump && lambda_fl > 0
        eta_heat = 1 - lambda_pump / lambda_fl;
    else
        warning(['apply_quantum_defect is true but fluorescence_wavelength is NaN or ', ...
                 'not > pump wavelength. Using eta_heat = 1 (all absorbed power is heat).']);
    end
end

% -------------------------------------------------------------------------
% 4. Super-Gaussian radial profile and Beer-Lambert heating
% -------------------------------------------------------------------------
phi_r = exp(-2 * (r / w).^order);          % normalised radial shape

integrand = 2 * pi * r .* phi_r;
P_shape = trapz(r, integrand);
if ~isfinite(P_shape) || P_shape <= 0
    error('Invalid pump normalisation integral. Check pump_radius_m, crystal_diameter_m, and nr.');
end

I0  = P_eff / P_shape;                     % [W/m²] peak on-axis intensity
I_r = I0 * phi_r;                          % [W/m²] radial profile

I_z  = exp(-alpha * z);                    % [—] axial attenuation
I_rz = I_r * I_z;                          % [W/m²] 2-D intensity field (nr × nz)
Q    = eta_heat * alpha * I_rz;            % [W/m³] volumetric heat source

% Fraction of pump power absorbed over crystal length
P_absorbed_frac = 1 - exp(-alpha * L);

% -------------------------------------------------------------------------
% 5. Saturation intensity
%    I_sat = h*nu_laser / (sigma * tau)
%    Uses emission_cross_section and interpolated lifetime at the user-set
%    doping level. Requires fluorescence_wavelength.
% -------------------------------------------------------------------------
I_sat = NaN;
tau = material.lifetime_vs_doping.interp(settings.Nd_concentration_wt_pct);
sigma = material.emission_cross_section;   % [m²]
lambda_fl_for_sat = material.fluorescence_wavelength;
if ~isnan(sigma) && ~isnan(tau) && ~isnan(lambda_fl_for_sat)
    nu_laser = SPEED_C / lambda_fl_for_sat;    % [Hz]
    I_sat = PLANCK_H * nu_laser / (sigma * tau);   % [W/m²]
end

% -------------------------------------------------------------------------
% 6. B-integral estimate (informational)
%    B = (2*pi/lambda_pump) * n2_SI * I_peak * L
%    Unit conversion: n2 [1e-13 e.s.u.] -> n2_SI [m²/W]
%    Using: n2_SI ≈ n2_esu * 1e-13 * (4.19e-7 / n0)  [approximate, Hellwarth 1977]
%    Flagged as estimate; provide measured n2 in SI if precision is needed.
% -------------------------------------------------------------------------
B_integral = NaN;
n2_esu_val = material.nonlinear_n2_esu;   % [1e-13 e.s.u.], upper bound from table
n0 = material.refractive_index_1053nm;
if ~isnan(n2_esu_val) && ~isnan(n0) && n0 > 0
    % Named conversion constant from Hellwarth (1977) for CGS e.s.u. to SI [m²/W]
    N2_CGS_TO_SI = 4.19e-7;   % [m²/W per e.s.u.], approximate (Hellwarth 1977)
    n2_SI = n2_esu_val * 1e-13 * (N2_CGS_TO_SI / n0);  % [m²/W] CGS-to-SI conversion
    B_integral = (2 * pi / lambda_pump) * n2_SI * I0 * L;
end

% -------------------------------------------------------------------------
% Assemble output
% -------------------------------------------------------------------------
pump = struct();
pump.r             = r;
pump.z             = z;
pump.alpha_808     = alpha;            % [1/m]
pump.I_rz          = I_rz;            % [W/m²]
pump.Q             = Q;               % [W/m³]
pump.eta_heat      = eta_heat;        % [—] quantum-defect fraction
pump.T_fresnel     = T_fresnel;       % [—] combined face transmittance
pump.P_absorbed_frac = P_absorbed_frac; % [—]
pump.I_sat         = I_sat;           % [W/m²] saturation intensity
pump.B_integral    = B_integral;      % [rad] B-integral, approximate
pump.lambda_pump_m = lambda_pump;     % [m]
end
