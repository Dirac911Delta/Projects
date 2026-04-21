% main.m
% Time-dependent thermal simulation for end-pumped Nd-doped laser glasses.
% MATLAB R2023 compatible.
%
% All material parameters are drawn from the three manufacturer data tables.
% Missing table entries are NaN; any gap that affects simulation is flagged.
% Refer to materials.m for the complete annotated database.

clear; clc;

%% ======================= User controls ================================ %%

% --- Material selection ------------------------------------------------ %
settings.material_name = 'N31';
% Valid: 'N31','N41','N51','NAP2','NAP4','NF1','NF2','NSG2'

% --- Pump and geometry ------------------------------------------------- %
settings.pump_wavelength_m       = 808e-9;    % [m]  pump diode wavelength
settings.pump_power_W            = 50;        % [W]
settings.pump_radius_m           = 0.6e-3;    % [m]  super-Gaussian 1/e radius
settings.supergaussian_order     = 6;         % super-Gaussian order n
settings.crystal_diameter_m      = 10e-3;     % [m]  fixed by problem statement
settings.crystal_thickness_m     = 20e-3;     % [m]  user-defined

% --- Doping ------------------------------------------------------------ %
settings.Nd_concentration_wt_pct = 1.0;       % [wt%] active doping level
                                               % used for lifetime lookup and
                                               % proportional alpha_808 estimate

% --- 808 nm absorption (NOT in manufacturer tables) -------------------- %
% Priority order: (1) user value  (2) material.absorption_808  (3) estimate
settings.absorption_808_user             = NaN;  % [1/m] measured value; set this
settings.estimate_absorption_from_doping = true; % fall back to proportional estimate
settings.absorption_808_ref              = NaN;  % [1/m] reference alpha at ref doping
settings.Nd_concentration_ref_wt_pct    = 1.0;  % [wt%] reference doping for above

% --- Physical effects flags -------------------------------------------- %
% Quantum-defect heating: Q = (1 - lambda_pump/lambda_fl) * alpha * I
% Uses material.fluorescence_wavelength. Recommended: true.
settings.apply_quantum_defect = true;

% Fresnel reflection at crystal faces (uses n(1053 nm) as proxy for n(808 nm))
% For AR-coated crystals set to false (default).
settings.apply_fresnel_reflection_loss = false;

% Fallback generic thermal props when table value is NaN.
% These are generic Nd-phosphate glass estimates; flagged in output.
settings.use_realistic_thermal_props = true;

% --- Boundary conditions ----------------------------------------------- %
settings.heat_sink_temperature_K = 293.15;        % [K]  293.15 = 20°C
settings.bc_z0.type  = 'insulated';               % 'insulated' | 'fixed'
settings.bc_z0.value_K = settings.heat_sink_temperature_K;
settings.bc_zL.type  = 'insulated';
settings.bc_zL.value_K = settings.heat_sink_temperature_K;

% --- Time and grid ----------------------------------------------------- %
settings.t_end_s = 2.0;    % [s]  total simulation time
settings.dt_s    = 0.01;   % [s]  time step
settings.nr      = 81;     % radial nodes
settings.nz      = 101;    % axial nodes

% --- Lifetime override ------------------------------------------------- %
settings.lifetime_override_s = NaN;  % [s] override interpolated value if set

%% ==================== Material and model setup ======================== %%
material = materials(settings.material_name);

% --- Lifetime (from doping table or user override) --------------------- %
if ~isnan(settings.lifetime_override_s)
    material.lifetime_s = settings.lifetime_override_s;
else
    material.lifetime_s = material.lifetime_vs_doping.interp(settings.Nd_concentration_wt_pct);
end

% --- Thermal property resolution --------------------------------------- %
%  Order: (1) direct field already in struct from table
%         (2) no override needed – table values are already filled by materials.m
% Apply user overrides only if explicitly set.

% (Optional) hard user overrides – leave NaN to use tabulated values.
thermal_conductivity_user = NaN;   % [W/(m·K)]
specific_heat_user        = NaN;   % [J/(kg·K)]
density_user              = NaN;   % [kg/m³]

if ~isnan(thermal_conductivity_user), material.thermal_conductivity = thermal_conductivity_user; end
if ~isnan(specific_heat_user),        material.specific_heat        = specific_heat_user;        end
if ~isnan(density_user),              material.density              = density_user;              end

% If a property is still NaN after table lookup, apply generic fallbacks
if settings.use_realistic_thermal_props
    if isnan(material.thermal_conductivity)
        material.thermal_conductivity = 0.80;   % [W/(m·K)] generic Nd-phosphate
        warning('ASSUMPTION: thermal_conductivity not in table for %s. Using 0.80 W/(m·K).', ...
                material.name);
    end
    if isnan(material.specific_heat)
        material.specific_heat = 800;           % [J/(kg·K)]
        warning('ASSUMPTION: specific_heat not in table for %s. Using 800 J/(kg·K).', ...
                material.name);
    end
    if isnan(material.density)
        material.density = 2700;                % [kg/m³]
        warning('ASSUMPTION: density not in table for %s. Using 2700 kg/m³.', ...
                material.name);
    end
end

if any(isnan([material.thermal_conductivity, material.specific_heat, material.density]))
    error(['Thermal properties incomplete for %s. ', ...
           'Provide values manually or enable settings.use_realistic_thermal_props.'], ...
           material.name);
end

if settings.crystal_thickness_m <= 0
    error('settings.crystal_thickness_m must be a positive value.');
end

% --- Validate absorption estimation inputs before calling pump_model --- %
if settings.estimate_absorption_from_doping && ...
   isnan(settings.absorption_808_user) && isnan(material.absorption_808)
    if any(isnan([settings.absorption_808_ref, ...
                  settings.Nd_concentration_ref_wt_pct, ...
                  settings.Nd_concentration_wt_pct]))
        error(['Cannot estimate alpha_808 for %s: set settings.absorption_808_ref, ', ...
               'Nd_concentration_ref_wt_pct, and Nd_concentration_wt_pct, ', ...
               'or provide settings.absorption_808_user.'], material.name);
    end
end

%% ==================== Run simulation ================================== %%
pump    = pump_model(material, settings);
results = heat_solver(material, pump, settings);
post    = post_processing(material, settings, results);

%% ==================== Console summary ================================= %%
fprintf('\n');
fprintf('=== Simulation Summary: %s ===\n', material.name);
fprintf('\n-- Material properties used --\n');
fprintf('  n(1053 nm)                    : %.4f\n',  material.refractive_index_1053nm);
fprintf('  dn/dT                         : %.4g 1/K\n', material.dn_dT);
fprintf('  Thermal OPL coeff (W0)        : %.4g 1/K\n', material.thermal_opl_coeff);
fprintf('  Thermal expansion coeff (CTE) : %.4g 1/K\n', material.thermal_expansion_coeff);
fprintf('  Thermal conductivity          : %.4f W/(m·K)\n', material.thermal_conductivity);
fprintf('  Specific heat                 : %.1f J/(kg·K)\n', material.specific_heat);
fprintf('  Density                       : %.0f kg/m³\n', material.density);
fprintf('  Fluorescence wavelength       : %.1f nm\n', material.fluorescence_wavelength * 1e9);
fprintf('  Emission cross section        : %.4g m²\n', material.emission_cross_section);
fprintf('  Effective bandwidth           : %.2f nm\n', material.effective_bandwidth * 1e9);
fprintf('  Lifetime at %.2f wt%%          : %.4g s\n', ...
        settings.Nd_concentration_wt_pct, material.lifetime_s);
fprintf('  Nd3+ concentration (tabulated): %.4g ions/m³\n', material.Nd_concentration);
fprintf('  Transition temperature        : %.0f °C\n', material.transition_temperature_C);
fprintf('  Softening temperature         : %.0f °C\n', material.softening_temperature_C);
fprintf('  Youngs modulus                : %.4g Pa\n', material.youngs_modulus_Pa);
fprintf('  Poisson''s ratio               : %.2f\n', material.poissons_ratio);
fprintf('  Fracture toughness            : %.2f MPa·m^0.5\n', material.fracture_toughness_MPa_sqrtm);
fprintf('  Nonlinear index n2            : %.4g x1e-13 e.s.u.\n', material.nonlinear_n2_esu);
fprintf('  Absorption coeff @1053 nm max : %.4g m⁻¹\n', material.absorption_coeff_1053nm_max);

fprintf('\n-- Pump model --\n');
fprintf('  alpha_808 used                : %.4g 1/m\n',  pump.alpha_808);
fprintf('  Quantum-defect eta_heat       : %.4f\n',      pump.eta_heat);
fprintf('  Fresnel face transmittance    : %.4f\n',      pump.T_fresnel);
fprintf('  Absorbed pump fraction        : %.4f (%.1f %%)\n', ...
        pump.P_absorbed_frac, pump.P_absorbed_frac * 100);
if ~isnan(pump.I_sat)
    fprintf('  Saturation intensity (I_sat)  : %.4g W/m²\n', pump.I_sat);
end
if ~isnan(pump.B_integral)
    fprintf('  B-integral (approx.)          : %.4f rad\n', pump.B_integral);
end

fprintf('\n-- Thermal results (t = %.3f s) --\n', results.t(end));
fprintf('  Peak temperature              : %.2f K  (%.1f °C)\n', ...
        post.peak_temperature_K(end), post.peak_temperature_C);
fprintf('  Max |dT/dr| at mid-plane      : %.4g K/m\n', post.radial_gradient_max_K_per_m(end));
if ~isnan(post.end_face_bulging_m)
    fprintf('  End-face axial bulging (axis) : %.4g m  (%.4g nm)\n', ...
            post.end_face_bulging_m, post.end_face_bulging_m * 1e9);
end

fprintf('\n-- Thermal lens --\n');
fprintf('  OPL coeff used (W0)           : %.4g 1/K\n', post.thermal_opl_coeff_used);
if ~isnan(post.thermal_lens_focal_power_1pm)
    if post.thermal_lens_focal_power_1pm ~= 0
        f_mm = 1e3 / post.thermal_lens_focal_power_1pm;
        fprintf('  Thermal lens focal power      : %.4g 1/m  (f = %.1f mm)\n', ...
                post.thermal_lens_focal_power_1pm, f_mm);
    else
        fprintf('  Thermal lens focal power      : 0 (no thermal gradient)\n');
    end
else
    fprintf('  Thermal lens focal power      : N/A (W0 and dn/dT both NaN)\n');
end

fprintf('\n-- Thermal stress (z_mid, t_end) --\n');
if ~isnan(post.max_tensile_stress_MPa)
    fprintf('  Max tensile stress (sigma_th) : %.4g MPa\n', post.max_tensile_stress_MPa);
    if ~isnan(material.fracture_toughness_MPa_sqrtm)
        fprintf('  Fracture toughness (K_IC)     : %.2f MPa·m^0.5\n', ...
                material.fracture_toughness_MPa_sqrtm);
        fprintf('  Note: Fracture risk requires crack size (a): K_I = sigma * sqrt(pi*a)\n');
    end
else
    fprintf('  Stress N/A (E, nu, or CTE missing for %s)\n', material.name);
end

fprintf('\n-- Temperature safety --\n');
switch post.temp_warning_level
    case 'none'
        fprintf('  OK – peak %.1f °C well below T_transition = %.0f °C\n', ...
                post.peak_temperature_C, material.transition_temperature_C);
    case 'CAUTION'
        fprintf('  CAUTION: peak %.1f °C within 50 °C of T_transition = %.0f °C\n', ...
                post.peak_temperature_C, material.transition_temperature_C);
    case 'CRITICAL'
        fprintf('  CRITICAL: peak %.1f °C EXCEEDS T_transition = %.0f °C! Reduce pump power.\n', ...
                post.peak_temperature_C, material.transition_temperature_C);
end

fprintf('\n-- Laser metrics --\n');
if ~isnan(post.small_signal_gain_dB)
    fprintf('  Small-signal gain G0          : %.2f dB  (g0·L = %.4f)\n', ...
            post.small_signal_gain_dB, post.gain_coeff_1pm * settings.crystal_thickness_m);
else
    fprintf('  Small-signal gain             : N/A (sigma_e or Nd_conc missing)\n');
end
if ~isnan(post.internal_transmission_1053nm)
    fprintf('  Internal transmission @1053nm : %.6f  (%.4f dB/m loss)\n', ...
            post.internal_transmission_1053nm, ...
            -10*log10(post.internal_transmission_1053nm)/settings.crystal_thickness_m);
end

fprintf('\n');
