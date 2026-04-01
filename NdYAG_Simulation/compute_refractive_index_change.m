function [delta_n, phi_TL, f_thermal] = compute_refractive_index_change(p, T_3d)
% COMPUTE_REFRACTIVE_INDEX_CHANGE  Compute the thermally and stress-optically
%   induced refractive index perturbation Δn(x,y,z) for [1 1 1]-cut Nd:YAG.
%
% Following the analytical model of Salinas-Alvarado et al. (ref [17] in
% Harrison 2023) and the parameters of Table 1, the total refractive-index
% change has two contributions:
%
%   Δn(r,z) = Δn_T(r,z)  +  Δn_SO(r,z)
%
% 1) Thermo-optic contribution:
%       Δn_T = (dn/dT) · ΔT(r,z)
%
% 2) Stress-optic (elasto-optic) contribution via thermoelastic stresses.
%    For an isotropic [1 1 1]-cut rod under plane-stress (free end-faces,
%    σ_z = 0) the thermoelastic stress state in (r,φ) is (Timoshenko &
%    Goodier, "Theory of Elasticity"):
%
%       σ_r(r)  = α_th·E/(1−ν) · [ ⟨ΔT⟩_R − ⟨ΔT⟩_r ]
%       σ_φ(r)  = α_th·E/(1−ν) · [ ⟨ΔT⟩_R + ⟨ΔT⟩_r − ΔT(r) ]
%
%    where ⟨ΔT⟩_R = (2/R²) ∫₀^R ΔT·r dr   (average over full cross-section)
%          ⟨ΔT⟩_r = (2/r²) ∫₀^r ΔT·r dr   (average within radius r)
%
%    Strains (plane stress):
%       ε_r = (σ_r − ν·σ_φ)/E + α_th·ΔT
%       ε_φ = (σ_φ − ν·σ_r)/E + α_th·ΔT
%       ε_z = −ν·(σ_r + σ_φ)/E + α_th·ΔT
%
%    For [1 1 1]-cut YAG, both polarisations experience the same Δn
%    (isotropic in the (111) plane), and the effective elasto-optic
%    coefficient is (Salinas-Alvarado 2012, J.Appl.Phys.111, 013112):
%       Δn_SO = −n₀³/2 · B_eff · (ε_r + ε_φ + ε_z)
%
%    where  B_eff = (p₁₁ + 2·p₁₂ + 4·p₄₄) / 3   for [111] propagation.
%
%    With the plane-stress thermoelastic result:
%       ε_r + ε_φ + ε_z = α_th·ΔT·(1 − 2ν) + (1−2ν)(σ_r + σ_φ)/E
%    (this is the volumetric dilation, θ = Tr(ε))
%
% 3) The accumulated phase shift along the full crystal length is then:
%       φ_TL(x,y) = k₀ · ∫₀^L Δn(x,y,z) dz
%    and the thermal-lens focal length (parabolic-phase approximation):
%       1/f = −(∂²φ_TL/∂r²)|_{r=0} / k₀
%           = −(∂²Δn_integrated/∂r²)|_{r=0}
%
% INPUTS
%   p       - parameter struct (NdYAG_params)
%   T_3d    - [N×N×Nz] temperature distribution           [K]
%
% OUTPUTS
%   delta_n - [N×N×Nz] refractive-index perturbation Δn    [ ]
%   phi_TL  - [N×N]    accumulated thermal-lens phase       [rad]
%   f_thermal - scalar thermal lens focal length            [m]
%              (positive → converging, negative → diverging)

N    = p.N_grid;
Nz   = p.N_z;
dz   = p.dz;
dx   = p.dx;
n0   = p.n0;
dndT = p.dndT;
ath  = p.alpha_th;
nu   = p.nu_poisson;
E_Y  = p.E_young;
p11  = p.p11;
p12  = p.p12;
p44  = p.p44;
k0s  = p.k_s;             % signal wavenumber in medium
T_s  = p.T_sink;
R_c  = p.R_c;
R_g  = p.R_grid;          % [N×N] radial coordinate map

% Effective elasto-optic coefficient for [111]-cut propagation:
%   B_eff = (p11 + 2·p12 + 4·p44) / 3
B_eff = (p11 + 2*p12 + 4*p44) / 3;

% ---- Allocate output ----------------------------------------------------
delta_n = zeros(N, N, Nz, 'single');

% ---- 1-D radial arrays for the thermoelastic integrals -----------------
Nr_1d   = N;
r_1d    = linspace(0, p.L_grid/2, Nr_1d)';  % [m]  column vector (N×1)
dr_1d   = r_1d(2) - r_1d(1);

% ---- Process each z-slice -----------------------------------------------
for kz = 1 : Nz
    DT = double(T_3d(:,:,kz)) - T_s;         % ΔT [K] [N×N]

    % -- Thermo-optic contribution ----------------------------------------
    Dn_T = dndT .* DT;

    % -- Stress-optic contribution (plane-stress thermoelastic) -----------
    % Work in cylindrical coordinates by azimuthal averaging of ΔT(r)
    % (valid since the initial fields are radially symmetric; BPM breaks
    %  symmetry only weakly for paraxial propagation)

    % Azimuthally averaged ΔT profile → ΔT_1d(r)
    DT_1d = zeros(Nr_1d, 1);
    for ir = 1:Nr_1d
        mask = (R_g >= r_1d(ir) - dr_1d/2) & (R_g < r_1d(ir) + dr_1d/2);
        if any(mask(:))
            DT_1d(ir) = mean(DT(mask));
        end
    end

    % Cumulative radial average ⟨ΔT⟩_r(r) = (2/r²) ∫₀^r ΔT·r' dr'
    % Avoid division by zero at r=0
    r_safe  = max(r_1d, dr_1d/4);
    int_rDT = cumtrapz(r_1d, DT_1d .* r_1d);        % ∫₀^r ΔT·r' dr'  [K·m²]
    DT_avg_r= 2 ./ r_safe.^2 .* int_rDT;             % ⟨ΔT⟩_r  [K]

    % Full cross-section average ⟨ΔT⟩_R (constant with r)
    % Integrate only up to R_c
    in1d     = r_1d <= R_c;
    DT_avg_R = 2/R_c^2 * trapz(r_1d(in1d), DT_1d(in1d) .* r_1d(in1d));

    % Thermoelastic stresses (MPa → Pa already in E_Y [Pa])
    prefactor = ath * E_Y / (1 - nu);
    sig_r_1d  = prefactor .* (DT_avg_R - DT_avg_r);    % [Pa]
    sig_phi_1d= prefactor .* (DT_avg_R + DT_avg_r - DT_1d);

    % Thermoelastic strains (plane stress)
    eps_r   = (sig_r_1d   - nu * sig_phi_1d) / E_Y + ath .* DT_1d;
    eps_phi = (sig_phi_1d - nu * sig_r_1d  ) / E_Y + ath .* DT_1d;
    eps_z   = -nu * (sig_r_1d + sig_phi_1d) / E_Y + ath .* DT_1d;

    % Volumetric dilation θ = ε_r + ε_φ + ε_z
    theta_1d = eps_r + eps_phi + eps_z;

    % Stress-optic index change on 1-D grid
    Dn_SO_1d = -n0^3 / 2 * B_eff .* theta_1d;

    % Interpolate 1-D Dn_SO back to 2-D grid
    Dn_SO = interp1(r_1d, Dn_SO_1d, R_g(:), 'linear', 0);
    Dn_SO = reshape(Dn_SO, N, N);

    % -- Total refractive index change ------------------------------------
    Dn_total = Dn_T + Dn_SO;
    Dn_total(R_g > R_c) = 0;       % only inside crystal

    delta_n(:,:,kz) = single(Dn_total);
end

% ---- Accumulated thermal-lens phase φ_TL(x,y) = k₀ ∫ Δn dz -----------
phi_TL = k0s * squeeze(sum(double(delta_n), 3)) * dz;    % [rad]  [N×N]

% ---- Thermal lens focal length (paraxial / parabolic approximation) ----
% Build a 1-D radial profile by azimuthally averaging phi_TL.
% Use the pre-computed |r| coordinate in p.R_grid (centred grid) to
% build a uniform radial sampling from 0 to R_c, independent of column
% ordering, to avoid issues with R_g(:,1) not passing through the origin.
Nr_tl = 256;
r_TL  = linspace(0, p.R_c * 1.2, Nr_tl)';   % uniform 1-D radial axis [m]
dr_tl = r_TL(2) - r_TL(1);
phi_r = zeros(Nr_tl, 1);
for ir = 1 : Nr_tl
    annulus = (R_g >= max(r_TL(ir) - dr_tl/2, 0)) & ...
              (R_g <  r_TL(ir) + dr_tl/2);
    if any(annulus(:))
        phi_r(ir) = mean(phi_TL(annulus));
    end
end
% Fit φ(r) = a + b·r² to the central region to get the parabolic coefficient
r_fit_max = min(p.R_c / 3, 5 * dx);    % fit within inner 1/3 of crystal
mask_fit  = r_TL <= r_fit_max;
if sum(mask_fit) >= 3
    r_fit  = r_TL(mask_fit);
    ph_fit = phi_r(mask_fit);
    coeffs = polyfit(r_fit.^2, ph_fit, 1);   % phi = coeffs(1)*r^2 + coeffs(2)
    b_par  = coeffs(1);   % d²φ/(2·dr²) at r=0 → d²φ/dr² = 2·b_par
    % 1/f_TL = −d²φ/dr² / k₀_vac = −2·b_par / k₀_vac
    k0_vac     = 2*pi / p.lambda_s;
    f_thermal  = -k0_vac / (2 * b_par);   % [m]
else
    f_thermal = Inf;
end

end
