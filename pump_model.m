function pump = pump_model(material, settings)
% pump_model.m
% End-pumped 808 nm model:
%   I(r) = I0 * exp(-2*(r/w)^6)
%   I(z) = I(r) * exp(-alpha*z)
%   Q(r,z) = alpha * I(r,z)
% where Q has units [W/m^3].

nr = settings.nr;
nz = settings.nz;
R = settings.crystal_diameter_m / 2;
L = settings.crystal_thickness_m;
w = settings.pump_radius_m;
P = settings.pump_power_W;
order = settings.supergaussian_order;

r = linspace(0, R, nr).';
z = linspace(0, L, nz);

% Determine alpha_808 (explicitly flagged missing-data path)
if ~isnan(settings.absorption_808_user)
    alpha = settings.absorption_808_user;
elseif ~isnan(material.absorption_808)
    alpha = material.absorption_808;
elseif settings.estimate_absorption_from_doping
    if isnan(settings.absorption_808_ref) || isnan(settings.Nd_concentration_ref_wt_pct) || isnan(settings.Nd_concentration_wt_pct)
        error('Cannot estimate alpha_808: missing reference or concentration inputs.');
    end
    % Assumption: alpha_808 proportional to Nd3+ concentration
    alpha = settings.absorption_808_ref * (settings.Nd_concentration_wt_pct / settings.Nd_concentration_ref_wt_pct);
    warning(['ASSUMPTION USED: alpha_808 estimated from doping proportionality. ', ...
             'Provide settings.absorption_808_user for measured value.']);
else
    error('alpha_808 unavailable. Set settings.absorption_808_user or enable estimate_absorption_from_doping.');
end

% Super-Gaussian radial shape
phi_r = exp(-2 * (r / w).^order);

% Normalize to pump power within simulation aperture
integrand = 2 * pi * r .* phi_r;
P_shape = trapz(r, integrand);
I0 = P / max(P_shape, eps);
I_r = I0 * phi_r;                         % [W/m^2]

% Axial Beer-Lambert attenuation and volumetric heating
I_z = exp(-alpha * z);                    % unitless
I_rz = I_r * I_z;                         % [W/m^2]
Q = alpha * I_rz;                         % [W/m^3]

pump = struct();
pump.r = r;
pump.z = z;
pump.alpha_808 = alpha;
pump.I_rz = I_rz;
pump.Q = Q;
end
