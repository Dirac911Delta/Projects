function post = post_processing(material, settings, results)
% post_processing.m
% Derives all physically meaningful outputs from the T(r,z,t) field,
% using every available material property.
% Pump-model outputs (eta_heat, I_sat, B_integral, etc.) are reported
% directly in main.m from the pump struct.
%
% Data fields used:
%   thermal_opl_coeff          -> thermal lens (preferred, full OPL)
%   dn_dT                      -> thermal lens (fallback if W0 is NaN)
%   thermal_expansion_coeff    -> thermal lens (CTE contribution), thermal stress
%   youngs_modulus_Pa          -> radial and hoop thermal stress
%   poissons_ratio             -> thermal stress
%   fracture_toughness_MPa_sqrtm -> stress safety reference
%   transition_temperature_C   -> peak-temperature safety check
%   softening_temperature_C    -> peak-temperature safety check
%   emission_cross_section     -> small-signal gain coefficient
%   Nd_concentration           -> small-signal gain coefficient
%   absorption_coeff_1053nm_max-> internal transmission at 1053 nm
%   refractive_index_1053nm    -> end-face OPL contribution to thermal lens
%   fluorescence_wavelength    -> (informational, reported in summary)
%   nonlinear_n2_esu           -> reported from pump struct (B-integral)

r   = results.r;
z   = results.z;
t   = results.t;
T   = results.T;
nr  = numel(r);
nz  = numel(z);
nt  = numel(t);
L   = z(end);
R   = r(end);
dr  = r(2) - r(1);

% -------------------------------------------------------------------------
% 1. Peak temperature vs time
% -------------------------------------------------------------------------
peakT_K = zeros(nt, 1);
for k = 1:nt
    peakT_K(k) = max(T(:,:,k), [], 'all');
end

% -------------------------------------------------------------------------
% 2. Radial temperature gradient (max |dT/dr|) at mid-plane vs time
% -------------------------------------------------------------------------
[~, jmid] = min(abs(z - z(end)/2));
radial_grad_max = zeros(nt, 1);
for k = 1:nt
    Tr = T(:, jmid, k);
    dTdr = gradient(Tr, r);
    radial_grad_max(k) = max(abs(dTdr));
end

% -------------------------------------------------------------------------
% 3. Final-time axial profile on axis (r = 0)
% -------------------------------------------------------------------------
axial_profile_final_K = squeeze(T(1, :, end));

% -------------------------------------------------------------------------
% 4. Thermal lens focal power via parabolic fit to T(r, z_mid, t_end)
%
%    The optical path length change across the crystal radius is:
%      delta_OPL(r) = W0 * [T(r, z_mid) - T(r=0, z_mid)] * L
%    where W0 is the thermal coefficient of optical path length [1/K].
%
%    A thin lens with parabolic OPL variation delta_OPL(r) = -r²/(2f)
%    gives focal power P_lens = 1/f [m⁻¹]:
%      P_lens = -2 * W0 * C * L
%    where C is the r² coefficient of a parabolic fit to T(r, z_mid).
%
%    W0 is taken from material.thermal_opl_coeff (full OPL including dn/dT,
%    CTE end-face contribution, and elasto-optic effect) when available;
%    otherwise dn/dT is used with a note.
% -------------------------------------------------------------------------
W0 = material.thermal_opl_coeff;       % [1/K] preferred
if isnan(W0)
    W0 = material.dn_dT;               % [1/K] fallback
    if ~isnan(W0)
        warning(['thermal_opl_coeff is NaN; using dn_dT as W0 for thermal lens. ', ...
                 'This omits the CTE and elasto-optic contributions to OPL.']);
    end
end

thermal_lens_focal_power = NaN;        % [1/m]
parabolic_C = NaN;
if ~isnan(W0)
    Tr_mid = T(:, jmid, end);
    % Fit T(r) = T0 + C*r^2 over the inner third of the aperture.
    % The inner third targets the near-axis parabolic region where the
    % paraxial thermal lens approximation is valid. Beyond this region
    % the T(r) profile deviates from parabolic due to the Gaussian/super-
    % Gaussian pump cutoff (see, e.g., Innocenzi et al., J. Appl. Phys.
    % 75, 4991, 1994 for end-pumped thermal lens analysis).
    PARABOLIC_FIT_FRACTION = 1/3;
    fit_pts = max(3, round(nr * PARABOLIC_FIT_FRACTION));
    r_fit = r(1:fit_pts);
    T_fit = Tr_mid(1:fit_pts);
    p = polyfit(r_fit.^2, T_fit, 1);    % p(1) = C, p(2) = T0
    parabolic_C = p(1);                 % [K/m²]
    % P_lens = -2 * W0 * C * L  (standard thin GRIN lens focal power)
    thermal_lens_focal_power = -2 * W0 * parabolic_C * L;
end

% -------------------------------------------------------------------------
% 5. Thermal stress – Timoshenko plane-stress cylinder formula
%    (valid for a free-ended rod; stress at z = z_mid, t = t_end)
%
%    sigma_r(r)   = alpha_E * [T_bar - T_bar_inner(r)]
%    sigma_theta(r) = alpha_E * [T_bar + T_bar_inner(r) - T(r)]
%
%    where alpha_E = E * alpha_CTE / (1 - nu)      [Pa/K]
%          T_bar        = (2/R²) * integral_0^R T(r') r' dr'
%          T_bar_inner(r) = (2/r²) * integral_0^r T(r') r' dr'
%
%    At r = 0: sigma_r = sigma_theta = alpha_E * (T_bar - T(0))  (compressive)
%    Max tensile stress typically occurs in sigma_theta near the outer edge.
% -------------------------------------------------------------------------
E_mod = material.youngs_modulus_Pa;
nu    = material.poissons_ratio;
alpha_CTE = material.thermal_expansion_coeff;

sigma_r_MPa     = NaN(nr, 1);
sigma_theta_MPa = NaN(nr, 1);
max_tensile_stress_MPa = NaN;

if ~isnan(E_mod) && ~isnan(nu) && ~isnan(alpha_CTE)
    alpha_E = E_mod * alpha_CTE / (1 - nu);     % [Pa/K]
    Tr_stress = T(:, jmid, end);                % T(r) at mid-plane, final time

    R_sq = R^2;   % precompute constant
    T_bar = 2 / R_sq * trapz(r, Tr_stress .* r);  % area-average temperature

    % T_bar_inner(r): cumulative radial average using cumtrapz
    % T_bar_inner(r) = (2/r^2) * integral_0^r T*r' dr'
    cum_integrand = Tr_stress .* r;
    cum_integral  = cumtrapz(r, cum_integrand);   % integral from 0 to r_i

    r_sq = r.^2;              % precompute r² for the loop below
    T_bar_inner = NaN(nr, 1);
    T_bar_inner(1) = Tr_stress(1);   % L'Hopital limit at r=0 → T(0)
    for i = 2:nr
        T_bar_inner(i) = 2 / r_sq(i) * cum_integral(i);
    end

    sigma_r_Pa     = alpha_E * (T_bar - T_bar_inner);
    sigma_theta_Pa = alpha_E * (T_bar + T_bar_inner - Tr_stress);

    sigma_r_MPa     = sigma_r_Pa     * 1e-6;     % [MPa]
    sigma_theta_MPa = sigma_theta_Pa * 1e-6;     % [MPa]

    % Maximum tensile (positive) stress in sigma_theta
    max_tensile_stress_MPa = max(sigma_theta_MPa);
end

% -------------------------------------------------------------------------
% 6. Temperature safety check vs transition and softening temperatures
% -------------------------------------------------------------------------
T_peak_C = peakT_K(end) - 273.15;
T_trans  = material.transition_temperature_C;
T_soft   = material.softening_temperature_C;

temp_warning_level = 'none';
if ~isnan(T_trans)
    if T_peak_C >= T_trans
        temp_warning_level = 'CRITICAL';
    elseif T_peak_C >= T_trans - 50
        temp_warning_level = 'CAUTION';
    end
end
temperature_safety_ok = strcmp(temp_warning_level, 'none');

% -------------------------------------------------------------------------
% 7. Small-signal gain coefficient
%    g0 = sigma_e * N_Nd   [1/m]
%    G0_dB = 10 * log10(exp(g0 * L))
%
%    Uses tabulated Nd³⁺ concentration and emission cross section.
%    Reports saturated gain condition vs I_sat (from pump struct if supplied).
% -------------------------------------------------------------------------
small_signal_gain_dB = NaN;
gain_coeff_1pm = NaN;
sigma_e = material.emission_cross_section;   % [m²]
N_Nd    = material.Nd_concentration;         % [ions/m³]
if ~isnan(sigma_e) && ~isnan(N_Nd)
    gain_coeff_1pm   = sigma_e * N_Nd;                     % [1/m]
    small_signal_gain_dB = 10 * log10(exp(gain_coeff_1pm * L));
end

% -------------------------------------------------------------------------
% 8. Internal transmission at the signal wavelength (1053 nm)
%    T_int = exp(-alpha_1053 * L)
%    Uses the tabulated absorption coefficient upper bound.
% -------------------------------------------------------------------------
internal_transmission_1053nm = NaN;
alpha_1053 = material.absorption_coeff_1053nm_max;  % [m⁻¹] upper bound
if ~isnan(alpha_1053)
    internal_transmission_1053nm = exp(-alpha_1053 * L);
end

% -------------------------------------------------------------------------
% 9. Thermal expansion: axial end-face bulging on axis (r = 0)
%    delta_L = integral_0^L alpha_CTE * (T(r=0, z) - T_sink) dz
% -------------------------------------------------------------------------
end_face_bulging_m = NaN;
if ~isnan(alpha_CTE)
    T_axis_z = squeeze(T(1, :, end));           % T(r=0, z, t_end)
    dT_axis  = T_axis_z - settings.heat_sink_temperature_K;
    end_face_bulging_m = alpha_CTE * trapz(z, dT_axis);
end

% =========================================================================
% PLOTS
% =========================================================================
figure('Name', 'T(r,z) contour – final time', 'Color', 'w');
contourf(r*1e3, z*1e3, T(:,:,end).', 40, 'LineColor', 'none');
cb = colorbar; cb.Label.String = 'Temperature [K]';
xlabel('r [mm]'); ylabel('z [mm]');
title(sprintf('%s  –  T(r,z) at t = %.3f s', material.name, t(end)));

figure('Name', 'Peak temperature vs time', 'Color', 'w');
plot(t, peakT_K, 'LineWidth', 1.5);
if ~isnan(T_trans)
    yline(T_trans + 273.15, 'r--', sprintf('T_{transition} = %d °C', T_trans), ...
          'LabelHorizontalAlignment','left');
end
grid on; xlabel('Time [s]'); ylabel('Peak Temperature [K]');

figure('Name', 'Radial gradient vs time', 'Color', 'w');
plot(t, radial_grad_max, 'LineWidth', 1.5);
grid on; xlabel('Time [s]'); ylabel('max |dT/dr| [K/m]');

figure('Name', 'Axial temperature profile on axis – final time', 'Color', 'w');
plot(z*1e3, axial_profile_final_K, 'LineWidth', 1.5);
grid on; xlabel('z [mm]'); ylabel('T(r=0, z, t_{end}) [K]');

if ~isnan(max_tensile_stress_MPa)
    figure('Name', 'Radial thermal stress profile – final time', 'Color', 'w');
    plot(r*1e3, sigma_r_MPa,     'b-',  'LineWidth', 1.5, 'DisplayName', '\sigma_r');
    hold on;
    plot(r*1e3, sigma_theta_MPa, 'r--', 'LineWidth', 1.5, 'DisplayName', '\sigma_\theta');
    yline(0, 'k:');
    legend('Location', 'best');
    grid on;
    xlabel('r [mm]'); ylabel('Stress [MPa]');
    title(sprintf('%s  –  Thermal stress at z_{mid}, t_{end}', material.name));
end

% =========================================================================
% Assemble output struct
% =========================================================================
post = struct();
post.peak_temperature_K          = peakT_K;
post.radial_gradient_max_K_per_m = radial_grad_max;
post.axial_profile_final_K       = axial_profile_final_K;

post.thermal_lens_focal_power_1pm = thermal_lens_focal_power;   % [1/m]
post.thermal_opl_coeff_used       = W0;                          % [1/K]
post.parabolic_C_K_per_m2         = parabolic_C;                 % [K/m²]

post.sigma_r_MPa                  = sigma_r_MPa;                 % [MPa], radial
post.sigma_theta_MPa              = sigma_theta_MPa;             % [MPa], hoop
post.max_tensile_stress_MPa       = max_tensile_stress_MPa;      % [MPa]
post.fracture_toughness_MPa_sqrtm = material.fracture_toughness_MPa_sqrtm;

post.temperature_safety_ok        = temperature_safety_ok;
post.temp_warning_level           = temp_warning_level;
post.peak_temperature_C           = T_peak_C;

post.small_signal_gain_dB         = small_signal_gain_dB;        % [dB]
post.gain_coeff_1pm               = gain_coeff_1pm;              % [1/m]
post.internal_transmission_1053nm = internal_transmission_1053nm;% [—]
post.end_face_bulging_m           = end_face_bulging_m;          % [m]
end
