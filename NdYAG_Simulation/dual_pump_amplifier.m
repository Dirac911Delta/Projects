function results = dual_pump_amplifier(P_pump_vec, P_seed, crystal_cfg, verbose)
% DUAL_PUMP_AMPLIFIER  Simulate a dual-end-pumped Nd:YAG rod amplifier.
%
% Extends the Harrison 2023 model to a symmetric dual-pump configuration
% where equal pump power enters both end-faces of the crystal simultaneously.
% This distributes the heat load and population inversion more uniformly
% along the crystal axis, which is especially beneficial for the user's
% longer crystal (L=70 mm, R=3 mm).
%
% CONFIGURATION
%   • Pump 1 (P_p/2): enters at face A (z=0), DOE-shaped FT→Gaussian
%   • Pump 2 (P_p/2): enters at face B (z=L), DOE-shaped FT→Gaussian
%   • Seed:  enters at face B (z=L), propagates contra-prop to pump 1
%
% PHYSICS
%   The total pump rate at each point is the incoherent sum of both pumps
%   (the two pump beams are from separate diodes → no interference):
%       Rp_total(r,z) = Rp_1(r,z) + Rp_2(r,z)
%       N₂(r,z) = N_tot · Rp_total / (Rp_total + 1/τ + σ_em·φ_s)
%   The heat source is similarly additive:
%       Q(r,z) = η_h · α_eff · [I_p1(r,z) + I_p2(r,z)]
%
% INPUTS
%   P_pump_vec  - [1×M] TOTAL pump power (P_p1 + P_p2, split 50:50) [W]
%                 default: [0 5 10 20 30 40 50 75 100]
%   P_seed      - seed power (default 0.1 W)                         [W]
%   crystal_cfg - 'user' (default, L=70 mm) or 'paper'
%   verbose     - print progress (default: true)
%
% OUTPUT  (same structure as single_pump_amplifier)

if nargin < 1 || isempty(P_pump_vec), P_pump_vec = [0 5 10 20 30 40 50 75 100]; end
if nargin < 2 || isempty(P_seed),     P_seed = [];        end
if nargin < 3 || isempty(crystal_cfg),crystal_cfg='user'; end
if nargin < 4 || isempty(verbose),    verbose = true;     end

p = NdYAG_params(crystal_cfg);
if ~isempty(P_seed), p.P_seed = P_seed; end

N  = p.N_grid;
Nz = p.N_z;
dz = p.dz;

alpha_eff = compute_absorption_xsec(p, 0);

if verbose
    fprintf('\n=== Dual-Pump Amplifier  [%s] ===\n', crystal_cfg);
    fprintf('Crystal: L=%.0f mm, R=%.0f mm, 0.5%% Nd:YAG\n', p.L_c*1e3, p.R_c*1e3);
    fprintf('alpha_eff = %.2f cm⁻¹\n', alpha_eff*0.01);
end

M         = numel(P_pump_vec);
P_out_vec = zeros(1,M);
gain_dB   = zeros(1,M);
f_TL_vec  = zeros(1,M);
T_max_vec = zeros(1,M);
I_p_last  = []; N2_last = []; T_last = [];

E_seed_B = make_gaussian_field(p, p.P_seed, p.omega_seed_B);

for m = 1:M
    Pp  = P_pump_vec(m);
    Pp1 = Pp / 2;   % pump 1 (face A)
    Pp2 = Pp / 2;   % pump 2 (face B)

    if verbose, fprintf('  P_total = %6.1f W (each side %.1f W) ... ', Pp, Pp1); tic; end

    % -- Pump 1: enters at z=0, propagates +z --------------------------
    E_p1 = compute_DOE_phase(p, Pp1, 'contra', 0);
    I_p1_3d = zeros(N,N,Nz,'single');
    for kz = 1:Nz
        I_p1_3d(:,:,kz) = single(abs(E_p1).^2);
        [E_p1,~] = bpm_propagate(E_p1, p, dz, p.lambda_p, p.n0, alpha_eff, 0, 0);
        E_p1 = E_p1 .* p.mask_crystal;
    end

    % -- Pump 2: enters at z=L, propagates −z --------------------------
    % Use the same DOE model but mirrored (pump enters from the other end)
    E_p2 = compute_DOE_phase(p, Pp2, 'contra', 0);
    I_p2_3d = zeros(N,N,Nz,'single');
    for kz = Nz:-1:1
        I_p2_3d(:,:,kz) = single(abs(E_p2).^2);
        [E_p2,~] = bpm_propagate(E_p2, p, dz, p.lambda_p, p.n0, alpha_eff, 0, 0);
        E_p2 = E_p2 .* p.mask_crystal;
    end

    % -- Total pump intensity (incoherent sum) -------------------------
    I_p_total_3d = I_p1_3d + I_p2_3d;   % [W/m²]

    % -- Heat source ---------------------------------------------------
    Q_vol = p.eta_h .* alpha_eff .* double(I_p_total_3d);

    % -- Thermal model ------------------------------------------------
    T_3d = solve_heat_2d(p, Q_vol, alpha_eff);

    % -- Refractive index change -------------------------------------
    [delta_n_3d, ~, f_TL] = compute_refractive_index_change(p, T_3d);

    % -- Population inversion (no signal first) ----------------------
    N2_3d = zeros(N,N,Nz,'single');
    for kz = 1:Nz
        [N2s,~] = compute_population_inversion(p, double(I_p_total_3d(:,:,kz)), zeros(N), alpha_eff);
        N2_3d(:,:,kz) = single(N2s);
    end

    % -- Propagate seed (contra to pump 1, i.e. B→A) -----------------
    for iter = 1:2
        E_s   = E_seed_B;
        I_s_3d = zeros(N,N,Nz,'single');
        for kz = Nz:-1:1
            I_s_3d(:,:,kz) = single(abs(E_s).^2);
            [~, g_s] = compute_population_inversion(p, ...
                double(I_p_total_3d(:,:,kz)), double(abs(E_s).^2), alpha_eff);
            dn_s = double(delta_n_3d(:,:,kz));
            [E_s,~] = bpm_propagate(E_s, p, -dz, p.lambda_s, p.n0, 0, dn_s, g_s);
            E_s = E_s .* p.mask_crystal;
        end
        for kz = 1:Nz
            [N2_3d(:,:,kz),~] = compute_population_inversion(p, ...
                double(I_p_total_3d(:,:,kz)), double(I_s_3d(:,:,kz)), alpha_eff);
        end
    end

    P_out       = sum(abs(E_s(:)).^2) * p.dx^2;
    P_out_vec(m)= P_out;
    if P_out > 0 && p.P_seed > 0
        gain_dB(m) = 10*log10(P_out / p.P_seed);
    end
    f_TL_vec(m)  = f_TL;
    T_max_vec(m) = max(T_3d(:));

    if m == M
        I_p_last = I_p_total_3d; N2_last = N2_3d; T_last = T_3d;
    end

    if verbose
        fprintf('P_out=%.3f W, G=%.1f dB, f_TL=%.0f mm, T_max=%.1f K  [%.1fs]\n',...
                P_out, gain_dB(m), f_TL*1e3, T_max_vec(m), toc);
    end
end

results.P_pump    = P_pump_vec;
results.P_out     = P_out_vec;
results.gain_dB   = gain_dB;
results.f_thermal = f_TL_vec;
results.T_max     = T_max_vec;
results.I_pump_3d = I_p_last;
results.N2_3d     = N2_last;
results.T_3d      = T_last;
results.p         = p;
results.alpha_eff = alpha_eff;

% ---- Comparison plot: dual vs single-end (qualitative) -----------------
plot_dual_results(results);

end

% =========================================================================
function E = make_gaussian_field(p, P, w0)
A0 = sqrt(2*P/(pi*w0^2));
E  = A0 .* exp(-p.R_grid.^2 ./ w0^2);
end

% =========================================================================
function draw_circle(r_mm, color)
theta = linspace(0, 2*pi, 360);
hold on;
plot(r_mm*cos(theta), r_mm*sin(theta), '--', 'Color', color, 'LineWidth', 1.2);
end

function plot_dual_results(res)
p = res.p;
figure('Name',sprintf('Dual-pump amplifier [%s]', p.label),...
       'NumberTitle','off','Position',[120 120 1200 800]);

subplot(2,3,1);
plot(res.P_pump, res.P_out,'b-o','LineWidth',1.5,'MarkerSize',6);
xlabel('Total pump power (W)'); ylabel('Output power (W)');
title('Output power (dual pump)'); grid on;

subplot(2,3,2);
plot(res.P_pump, res.gain_dB,'r-s','LineWidth',1.5,'MarkerSize',6);
xlabel('Total pump power (W)'); ylabel('Gain (dB)');
title('Amplifier gain'); grid on;

subplot(2,3,3);
valid = isfinite(res.f_thermal) & res.f_thermal ~= 0;
if any(valid)
    plot(res.P_pump(valid), res.f_thermal(valid)*1e3,'g-d','LineWidth',1.5,'MarkerSize',6);
    xlabel('Total pump power (W)'); ylabel('f_{thermal} (mm)');
    title('Thermal lens focal length'); grid on;
end

if ~isempty(res.I_pump_3d)
    x_mm = p.x_vec * 1e3;
    subplot(2,3,4);
    imagesc(x_mm, x_mm, double(res.I_pump_3d(:,:,1))/max(max(double(res.I_pump_3d(:,:,1)))));
    colorbar; colormap(gca,'hot'); axis image;
    title('Norm. total pump at z=0'); xlabel('x (mm)'); ylabel('y (mm)');
    draw_circle(p.R_c*1e3, 'w');

    subplot(2,3,5);
    % Axial profile of on-axis N₂
    cx = floor(p.N_grid/2) + 1;   % centre pixel index (works for even and odd N)
    N2_axis = squeeze(double(res.N2_3d(cx, cx, :)));
    plot(p.z_vec*1e3, N2_axis*1e-6,'m-','LineWidth',1.5);
    xlabel('z (mm)'); ylabel('N_2 (×10^{12} cm^{-3})');
    title('On-axis population inversion N_2(z)'); grid on;

    subplot(2,3,6);
    T_axis = squeeze(double(res.T_3d(cx, cx, :)));
    plot(p.z_vec*1e3, T_axis,'k-','LineWidth',1.5);
    xlabel('z (mm)'); ylabel('T (K)');
    title('On-axis temperature T(z)'); grid on;
end

sgtitle(sprintf('Dual-pump amplifier: %s', p.label),'FontWeight','bold');
drawnow;
end
