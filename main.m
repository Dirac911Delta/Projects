% main.m
% Time-dependent thermal simulation for end-pumped Nd-doped laser glasses.
% MATLAB R2023 compatible.
%
% NOTE ON DATA:
% - Missing values from source tables are represented as NaN.
% - 808 nm absorption is not provided in the available optical table; this
%   script allows user input or a clearly-flagged concentration-scaled estimate.

clear; clc;

%% --------------------------- User controls ---------------------------- %%
settings.material_name = 'N31';            % {'N31','N41','N51','NAP2','NAP4','NF1','NF2','NSG2'}
settings.use_realistic_thermal_props = true;

% Pump and geometry
settings.pump_power_W = 50;                % [W]
settings.pump_radius_m = 0.6e-3;           % super-Gaussian radius w [m]
settings.supergaussian_order = 6;          % required by statement
settings.crystal_diameter_m = 10e-3;       % fixed by statement [m]
settings.crystal_thickness_m = 20e-3;      % user-defined [m]

% Doping inputs (for lifetime lookup/override and absorption estimate)
settings.Nd_concentration_wt_pct = 1.0;    % user-set [wt%]

% 808 nm absorption handling (CRITICAL missing data path)
settings.absorption_808_user = NaN;        % [1/m], user override if known
settings.estimate_absorption_from_doping = true;
settings.absorption_808_ref = NaN;         % [1/m], REQUIRED if using proportional estimate
settings.Nd_concentration_ref_wt_pct = 1.0;% [wt%], reference for proportional estimate

% Thermal properties (if not in material table)
settings.thermal_conductivity_user = NaN;  % [W/(m*K)]
settings.specific_heat_user = NaN;         % [J/(kg*K)]
settings.density_user = NaN;               % [kg/m^3]

% Boundary conditions
settings.heat_sink_temperature_K = 293.15; % [K]
settings.bc_z0.type = 'insulated';         % {'insulated','fixed'}
settings.bc_z0.value_K = settings.heat_sink_temperature_K;
settings.bc_zL.type = 'insulated';         % {'insulated','fixed'}
settings.bc_zL.value_K = settings.heat_sink_temperature_K;

% Time/grid controls
settings.t_end_s = 2.0;                    % [s]
settings.dt_s = 0.01;                      % [s]
settings.nr = 81;                          % radial nodes
settings.nz = 101;                         % axial nodes

% Optional override for lifetime model
settings.lifetime_override_s = NaN;        % [s]

%% ---------------------- Material and model setup --------------------- %%
material = materials(settings.material_name);

% Lifetime override or lookup
if ~isnan(settings.lifetime_override_s)
    material.lifetime_s = settings.lifetime_override_s;
else
    material.lifetime_s = material.lifetime_vs_doping.interp(settings.Nd_concentration_wt_pct);
end

% Thermal property handling
if ~isnan(settings.thermal_conductivity_user), material.thermal_conductivity = settings.thermal_conductivity_user; end
if ~isnan(settings.specific_heat_user), material.specific_heat = settings.specific_heat_user; end
if ~isnan(settings.density_user), material.density = settings.density_user; end

if settings.use_realistic_thermal_props
    % Clearly flagged generic defaults for Nd-doped glass when not supplied
    if isnan(material.thermal_conductivity), material.thermal_conductivity = 0.80; end  % [W/(m*K)]
    if isnan(material.specific_heat), material.specific_heat = 800; end                 % [J/(kg*K)]
    if isnan(material.density), material.density = 2700; end                            % [kg/m^3]
end

if any(isnan([material.thermal_conductivity, material.specific_heat, material.density]))
    error('Thermal properties are incomplete. Provide user values or enable realistic defaults.');
end

if settings.estimate_absorption_from_doping && isnan(settings.absorption_808_user) && isnan(material.absorption_808)
    if any(isnan([settings.absorption_808_ref, settings.Nd_concentration_ref_wt_pct, settings.Nd_concentration_wt_pct]))
        error(['Missing absorption-estimation inputs in main.m: set absorption_808_ref, ', ...
               'Nd_concentration_ref_wt_pct, and Nd_concentration_wt_pct, or provide absorption_808_user.']);
    end
end

pump = pump_model(material, settings);
results = heat_solver(material, pump, settings);
post = post_processing(material, settings, results);

%% ------------------------------ Console ------------------------------ %%
fprintf('\nMaterial: %s\n', material.name);
fprintf('alpha_808 used: %.6g 1/m\n', pump.alpha_808);
fprintf('Final peak temperature: %.2f K\n', post.peak_temperature_K(end));
fprintf('Estimated thermal lens power proxy (dn/dT * dT/dr max): %.3e 1/m\n', post.thermal_lens_proxy_per_m);
