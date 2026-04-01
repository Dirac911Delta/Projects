function results = single_pump_amplifier(P_pump_vec, P_seed, prop_dir, crystal_cfg, verbose)
% SINGLE_PUMP_AMPLIFIER  Simulate a single-pass end-pumped Nd:YAG rod
%   amplifier using the 3-D analytical / iterative-Fourier-BPM model of
%   Harrison, Forbes & Naidoo, Opt. Express 31(15), 24516 (2023).
%
% This function reproduces the paper's baseline results (Figures 3-7) for
% the contra-propagating MOPA configuration:
%   • Pump enters at face A (z = 0) as a DOE-shaped beam.
%   • Seed enters at face B (z = L) and propagates in −z direction.
%   • Crystal: L = 25 mm, R = 2 mm, 0.5 at.% Nd:YAG (Table 1).
%
% ALGORITHM (Section 3.5, Harrison 2023):
%   For each pump power P_pump:
%   1.  Generate initial pump field Ψ_G at z=0 using the DOE phase
%       (Equations 1-4).
%   2.  Propagate the pump in +z steps via the angular-spectrum BPM,
%       absorbing power at each slice (Beer-Lambert with α_eff).
%   3.  Store local pump intensity I_p(x,y,z) at every slice.
%   4.  Solve the 3-D heat equation to get T(x,y,z) and hence Δn(x,y,z).
%   5.  Compute steady-state population inversion N₂(x,y,z) from the
%       4-level rate equations.
%   6.  Generate the seed field at z = L and propagate it in −z steps
%       with amplitude gain g = σ_em · N₂ and thermal-lens phase from Δn.
%   7.  Record output seed power, beam quality, and thermal-lens focal
%       length.
%
% INPUTS
%   P_pump_vec  - [1×M] pump power values to sweep                    [W]
%                 (default: [0, 5, 10, 15, 20, 25, 30, 35, 38])
%   P_seed      - seed input power (default from params: 0.1 W)       [W]
%   prop_dir    - 'contra' (default, paper configuration) or 'co'
%   crystal_cfg - 'paper' (default, L=25 mm) or 'user' (L=70 mm)
%   verbose     - true/false – print progress (default: true)
%
% OUTPUT
%   results  struct with fields:
%     .P_pump        [W]   pump power vector
%     .P_out         [W]   output seed power after amplification
%     .gain_dB       [dB]  single-pass gain 10·log10(P_out/P_seed)
%     .f_thermal     [m]   thermal lens focal length vs P_pump
%     .T_max         [K]   peak crystal temperature vs P_pump
%     .I_pump_3d     {cell} I_p(x,y,z) for last pump power (for plotting)
%     .N2_3d         {cell} N₂(x,y,z) for last pump power
%     .T_3d          {cell} T(x,y,z) for last pump power
%     .p             parameter struct used

if nargin < 1 || isempty(P_pump_vec), P_pump_vec = [0 5 10 15 20 25 30 35 38]; end
if nargin < 2 || isempty(P_seed),     P_seed = [];          end
if nargin < 3 || isempty(prop_dir),   prop_dir = 'contra';  end
if nargin < 4 || isempty(crystal_cfg),crystal_cfg = 'paper';end
if nargin < 5 || isempty(verbose),    verbose = true;       end

% ---- Parameters ---------------------------------------------------------
p = NdYAG_params(crystal_cfg);
if ~isempty(P_seed), p.P_seed = P_seed; end

N   = p.N_grid;
Nz  = p.N_z;
dz  = p.dz;

% Effective absorption cross-section (spectral, no diode heating offset)
alpha_eff = compute_absorption_xsec(p, 0);

if verbose
    fprintf('\n=== Single-Pump Amplifier  [%s]  (%s) ===\n', crystal_cfg, prop_dir);
    fprintf('Crystal: L=%.0f mm, R=%.0f mm, 0.5%% Nd:YAG\n', ...
            p.L_c*1e3, p.R_c*1e3);
    fprintf('alpha_eff = %.2f cm⁻¹\n', alpha_eff*0.01);
    fprintf('I_sat     = %.2e W/m²\n', p.I_sat);
end

% ---- Allocate output arrays --------------------------------------------
M          = numel(P_pump_vec);
P_out_vec  = zeros(1, M);
gain_dB    = zeros(1, M);
f_TL_vec   = zeros(1, M);
T_max_vec  = zeros(1, M);

% Stores for the last pump power (for plotting)
I_pump_last = [];
N2_last     = [];
T_last      = [];

% ---- Seed field at crystal entry face ----------------------------------
% Gaussian TEM₀₀ seed, waist at z=0 (face A) in contra-prop or at z=L in co-prop
%   In contra-prop: seed enters at face B (z=L), waist at z=0 (face A)
%   → at z=L: seed radius = ω_seed_B (from Table 1)
E_seed_B = make_gaussian_field(p, p.P_seed, p.omega_seed_B);

% ---- Sweep over pump powers --------------------------------------------
for m = 1:M
    Pp = P_pump_vec(m);
    if verbose
        fprintf('  P_pump = %5.1f W ... ', Pp);
        tic;
    end

    % == Step 1: Initial pump field at crystal face A (z=0) ==============
    E_pump = compute_DOE_phase(p, Pp, prop_dir, 0);

    % == Step 2: Propagate pump through crystal (+z) =====================
    I_pump_3d = zeros(N, N, Nz, 'single');
    E_p       = E_pump;

    for kz = 1 : Nz
        I_pump_3d(:,:,kz) = single(abs(E_p).^2);

        % Half-step: record intensity BEFORE propagation
        % Propagate pump with absorption (no gain, no thermal phase for pump)
        [E_p, ~] = bpm_propagate(E_p, p, dz, p.lambda_p, p.n0, ...
                                  alpha_eff, 0, 0);

        % Clip to crystal aperture (hard wall)
        E_p = E_p .* p.mask_crystal;
    end

    % == Step 3: Compute heat source Q = eta_h · alpha_eff · I_p  =========
    Q_vol = p.eta_h .* alpha_eff .* double(I_pump_3d);   % [W/m³]

    % == Step 4: Solve heat equation for T(x,y,z) =========================
    T_3d = solve_heat_2d(p, Q_vol, alpha_eff);

    % == Step 5: Refractive-index change Δn(x,y,z) ========================
    [delta_n_3d, phi_TL, f_TL] = compute_refractive_index_change(p, T_3d);

    % == Step 6: Population inversion N₂(x,y,z) ===========================
    N2_3d = zeros(N, N, Nz, 'single');
    for kz = 1:Nz
        I_s_2d = zeros(N);  % first pass: no signal yet (small-signal)
        [N2_slice, ~] = compute_population_inversion(p, ...
            double(I_pump_3d(:,:,kz)), I_s_2d, alpha_eff);
        N2_3d(:,:,kz) = single(N2_slice);
    end

    % == Step 7: Propagate seed contra-propagating (B→A, i.e. z=L→0) ======
    E_s = E_seed_B;

    % Iterative scheme: two passes to account for gain saturation
    for iter = 1:2
        E_s = E_seed_B;  % restart from input each iteration
        I_s_3d = zeros(N, N, Nz, 'single');

        for kz = Nz : -1 : 1
            I_s_3d(:,:,kz) = single(abs(E_s).^2);

            % Gain at this slice (saturated by signal)
            [~, g_slice] = compute_population_inversion(p, ...
                double(I_pump_3d(:,:,kz)), double(abs(E_s).^2), alpha_eff);

            % Thermal phase per slice
            dn_slice = double(delta_n_3d(:,:,kz));

            % BPM step for backward (−z) propagation.
            % In the paraxial regime the gain profile and phase modulation
            % are applied locally at each slice; only the free-space
            % transfer function H = exp(i·kz·dz) changes sign for −z.
            % We negate dz so that H → exp(−i·kz·dz), correctly
            % implementing backward propagation through each slice.
            [E_s, ~] = bpm_propagate(E_s, p, -dz, p.lambda_s, p.n0, ...
                                      0, dn_slice, g_slice);
            E_s = E_s .* p.mask_crystal;
        end

        % Update N₂ with saturated signal for next iteration
        for kz = 1:Nz
            [N2_3d(:,:,kz), ~] = compute_population_inversion(p, ...
                double(I_pump_3d(:,:,kz)), double(I_s_3d(:,:,kz)), alpha_eff);
        end
    end

    % == Output quantities ================================================
    I_out = abs(E_s).^2;
    P_out = sum(I_out(:)) * p.dx^2;     % [W] output power at face A
    P_out_vec(m) = P_out;
    if P_out > 0 && p.P_seed > 0
        gain_dB(m) = 10 * log10(P_out / p.P_seed);
    end
    f_TL_vec(m)  = f_TL;
    T_max_vec(m) = max(T_3d(:));

    % Save full 3-D fields for the last pump power (plotting)
    if m == M
        I_pump_last = I_pump_3d;
        N2_last     = N2_3d;
        T_last      = T_3d;
    end

    if verbose
        fprintf('P_out = %.3f W, G = %.1f dB, f_TL = %.0f mm, T_max = %.1f K  [%.1fs]\n',...
                P_out, gain_dB(m), f_TL*1e3, T_max_vec(m), toc);
    end
end

% ---- Pack results -------------------------------------------------------
results.P_pump    = P_pump_vec;
results.P_out     = P_out_vec;
results.gain_dB   = gain_dB;
results.f_thermal = f_TL_vec;
results.T_max     = T_max_vec;
results.I_pump_3d = I_pump_last;
results.N2_3d     = N2_last;
results.T_3d      = T_last;
results.p         = p;
results.alpha_eff = alpha_eff;

% ---- Plot results -------------------------------------------------------
plot_amplifier_results(results, sprintf('Single-pump amplifier [%s, %s]', ...
                       crystal_cfg, prop_dir));

end

% =========================================================================
%  Local helper: draw crystal boundary circle on current axes
% =========================================================================
function draw_circle(r_mm, color)
theta = linspace(0, 2*pi, 360);
hold on;
plot(r_mm*cos(theta), r_mm*sin(theta), '--', 'Color', color, 'LineWidth', 1.2);
end

% =========================================================================
%  Local helper: make a Gaussian TEM₀₀ field
% =========================================================================
function E = make_gaussian_field(p, P, w0)
% Gaussian field amplitude normalised so that ∫|E|² dA = P
R  = p.R_grid;
A0 = sqrt(2*P / (pi * w0^2));
E  = A0 .* exp(-R.^2 ./ w0^2);
end

% =========================================================================
%  Local helper: plot amplifier results
% =========================================================================
function plot_amplifier_results(res, fig_title)
p  = res.p;
Pp = res.P_pump;
figure('Name', fig_title, 'NumberTitle','off', 'Position',[100 100 1400 900]);

% -- (a) Output power vs pump power --------------------------------------
subplot(2,3,1);
if p.P_seed > 0
    plot(Pp, res.P_out, 'b-o', 'LineWidth',1.5, 'MarkerSize',6);
    xlabel('Pump power (W)'); ylabel('Output power (W)');
    title('Output power vs pump power');
    grid on;
end

% -- (b) Single-pass gain ------------------------------------------------
subplot(2,3,2);
plot(Pp, res.gain_dB, 'r-s', 'LineWidth',1.5, 'MarkerSize',6);
xlabel('Pump power (W)'); ylabel('Gain (dB)');
title('Single-pass gain'); grid on;

% -- (c) Thermal lens focal length ---------------------------------------
subplot(2,3,3);
valid = isfinite(res.f_thermal) & (res.f_thermal ~= 0);
if any(valid)
    plot(Pp(valid), res.f_thermal(valid)*1e3, 'g-d','LineWidth',1.5,'MarkerSize',6);
    xlabel('Pump power (W)'); ylabel('f_{thermal} (mm)');
    title('Thermal lens focal length'); grid on;
end

% -- (d) Pump intensity at crystal input (z=0) ---------------------------
if ~isempty(res.I_pump_3d)
    subplot(2,3,4);
    I_in = double(res.I_pump_3d(:,:,1));
    x_mm = p.x_vec * 1e3;
    imagesc(x_mm, x_mm, I_in / max(I_in(:)));
    colorbar; colormap(gca, 'hot'); axis image;
    xlabel('x (mm)'); ylabel('y (mm)');
    title(sprintf('Norm. pump intensity at z=0\n(P_p=%.0f W)', res.P_pump(end)));
    draw_circle(p.R_c*1e3, 'w');

    % -- (e) Population inversion at central z-slice ---------------------
    subplot(2,3,5);
    kz_mid = round(p.N_z / 4);   % 1/4 into crystal (near pump face)
    N2_slice = double(res.N2_3d(:,:,kz_mid));
    imagesc(x_mm, x_mm, N2_slice * 1e-6);   % convert to cm⁻³ * 1e12
    colorbar; colormap(gca,'parula'); axis image;
    xlabel('x (mm)'); ylabel('y (mm)');
    title(sprintf('N_2 (×10^{12} cm^{-3}) at z=%.1f mm', p.z_vec(kz_mid)*1e3));
    draw_circle(p.R_c*1e3, 'w');

    % -- (f) Temperature distribution at central z-slice ----------------
    subplot(2,3,6);
    T_slice = double(res.T_3d(:,:,kz_mid));
    imagesc(x_mm, x_mm, T_slice);
    colorbar; colormap(gca,'jet'); axis image;
    xlabel('x (mm)'); ylabel('y (mm)');
    title(sprintf('Temperature (K) at z=%.1f mm', p.z_vec(kz_mid)*1e3));
    draw_circle(p.R_c*1e3, 'w');
end

sgtitle(fig_title, 'FontWeight','bold');
drawnow;
end
