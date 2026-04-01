function results = CW_flat_flat_laser(P_pump_vec, crystal_cfg, verbose)
% CW_FLAT_FLAT_LASER  Simulate a CW Nd:YAG laser in a flat/flat cavity.
%
% Extends the Harrison 2023 pump-mode-evolution model to a CW oscillator.
% A flat/flat (plano-plano) resonator is formed by placing a high-reflector
% (HR) at face A (z=0) and an output coupler (OC) at face B (z=L).
%
% PHYSICS
% -------
% The intracavity field builds up from spontaneous emission noise until the
% round-trip gain equals round-trip loss (threshold condition):
%
%   Round-trip gain   : G_rt = exp(2 · ∫₀^L g(r,z) dz)   [power]
%   Round-trip losses : L_rt = (1−R_OC) · (1−R_HR) · (1 + loss_other)
%   Threshold        : G_rt · (1−L_rt) = 1
%
% Above threshold, the steady-state intracavity power P_ic is found by
% iterating the Fox-Li round-trip propagation until convergence:
%
%   1. Propagate seed field through crystal (z=0→L) with gain and
%      thermal-lens phase using the BPM (same as amplifier model).
%   2. Apply OC reflectivity: E → √R_OC · E
%   3. Propagate back (z=L→0) with gain and thermal phase.
%   4. Apply HR reflectivity: E → √R_HR · E
%   5. Repeat until |P_ic(n+1)/P_ic(n) − 1| < tolerance.
%
% Output power:
%   P_out = (1 − R_OC) · P_ic
%
% Slope efficiency (above threshold, simplified):
%   η_slope = η_OC · η_overlap · η_pump_abs · η_Stokes
%
% INPUTS
%   P_pump_vec  - [1×M] pump power                                     [W]
%                 default: [0 2 5 8 10 15 20 25 30 35 38]
%   crystal_cfg - 'paper' (L=25 mm) or 'user' (L=70 mm) [default 'paper']
%   verbose     - print progress (default true)
%
% OUTPUT
%   results  struct with:
%     .P_pump          pump power vector                               [W]
%     .P_out           output power                                    [W]
%     .P_intracavity   intracavity power                               [W]
%     .P_threshold     estimated threshold pump power                  [W]
%     .slope_eff       slope efficiency above threshold                [W/W]
%     .f_thermal       thermal lens focal length                       [m]
%     .T_max           peak crystal temperature                        [K]
%     .p               parameter struct

if nargin < 1 || isempty(P_pump_vec), P_pump_vec = [0 2 5 8 10 15 20 25 30 35 38]; end
if nargin < 2 || isempty(crystal_cfg), crystal_cfg = 'paper'; end
if nargin < 3 || isempty(verbose),     verbose = true;        end

p = NdYAG_params(crystal_cfg);

N  = p.N_grid;
Nz = p.N_z;
dz = p.dz;

R_OC   = p.R_OC;
R_HR   = p.R_HR;
L_rt   = (1 - R_OC) * (1 - R_HR) + p.loss_rt_other;   % total round-trip loss

alpha_eff = compute_absorption_xsec(p, 0);

if verbose
    fprintf('\n=== CW Flat/Flat Laser  [%s] ===\n', crystal_cfg);
    fprintf('R_OC = %.2f, R_HR = %.4f, RT loss = %.3f\n', R_OC, R_HR, L_rt);
    fprintf('alpha_eff = %.2f cm⁻¹\n', alpha_eff*0.01);
end

M          = numel(P_pump_vec);
P_out_vec  = zeros(1,M);
P_ic_vec   = zeros(1,M);
f_TL_vec   = zeros(1,M);
T_max_vec  = zeros(1,M);

% Fox-Li convergence parameters
max_iter   = 150;
tol        = 1e-4;

% Initialise intracavity field as Gaussian seed (low power = noise level)
P_ic_init  = 1e-6;   % 1 µW noise seed to start Fox-Li

for m = 1:M
    Pp = P_pump_vec(m);
    if verbose, fprintf('  P_pump = %5.1f W ... ', Pp); tic; end

    % == Pump propagation (+z) ===========================================
    E_pump = compute_DOE_phase(p, Pp, 'contra', 0);
    I_p_3d = zeros(N,N,Nz,'single');
    E_p    = E_pump;
    for kz = 1:Nz
        I_p_3d(:,:,kz) = single(abs(E_p).^2);
        [E_p,~] = bpm_propagate(E_p, p, dz, p.lambda_p, p.n0, alpha_eff, 0, 0);
        E_p = E_p .* p.mask_crystal;
    end

    % == Thermal model ===================================================
    Q_vol    = p.eta_h .* alpha_eff .* double(I_p_3d);
    T_3d     = solve_heat_2d(p, Q_vol, alpha_eff);
    [delta_n_3d, ~, f_TL] = compute_refractive_index_change(p, T_3d);

    % == Small-signal gain check =========================================
    % Compute maximum on-axis gain integral to estimate threshold
    g_on_axis = zeros(Nz,1);
    for kz = 1:Nz
        [~,g_sl] = compute_population_inversion(p, ...
            double(I_p_3d(:,:,kz)), zeros(N), alpha_eff);
        cx = floor(N/2) + 1;   % centre pixel (valid for even and odd N)
        cy = cx;
        g_on_axis(kz) = g_sl(cy,cx);
    end
    G_ss_power = exp(2 * sum(g_on_axis) * dz);  % round-trip small-signal gain
    above_threshold = G_ss_power > (1 + L_rt);

    % == Fox-Li iteration (only above threshold) =========================
    P_ic = P_ic_init;

    if above_threshold && Pp > 0
        % Initialise intracavity field: Gaussian mode
        E_ic = make_gaussian_field(p, P_ic, p.omega_seed_A);

        for iter = 1:max_iter
            P_ic_prev = sum(abs(E_ic(:)).^2) * p.dx^2;

            % -- Forward pass A→B ----------------------------------------
            E_fwd = E_ic;
            for kz = 1:Nz
                [~,g_s] = compute_population_inversion(p, ...
                    double(I_p_3d(:,:,kz)), double(abs(E_fwd).^2), alpha_eff);
                dn_s = double(delta_n_3d(:,:,kz));
                [E_fwd,~] = bpm_propagate(E_fwd, p, dz, p.lambda_s, p.n0, 0, dn_s, g_s);
                E_fwd = E_fwd .* p.mask_crystal;
            end
            % Apply OC (partial transmission)
            E_fwd = E_fwd .* sqrt(R_OC);

            % -- Backward pass B→A (use −dz for correct backward propagation) --
            E_bwd = E_fwd;
            for kz = Nz:-1:1
                [~,g_s] = compute_population_inversion(p, ...
                    double(I_p_3d(:,:,kz)), double(abs(E_bwd).^2), alpha_eff);
                dn_s = double(delta_n_3d(:,:,kz));
                [E_bwd,~] = bpm_propagate(E_bwd, p, -dz, p.lambda_s, p.n0, 0, dn_s, g_s);
                E_bwd = E_bwd .* p.mask_crystal;
            end
            % Apply HR
            E_bwd = E_bwd .* sqrt(R_HR);
            % Apply additional round-trip loss
            E_bwd = E_bwd .* sqrt(1 - p.loss_rt_other);

            % -- Check convergence ----------------------------------------
            P_ic_new = sum(abs(E_bwd(:)).^2) * p.dx^2;
            if P_ic_prev > 0 && abs(P_ic_new/P_ic_prev - 1) < tol
                E_ic = E_bwd;
                break;
            end

            % -- Normalise to prevent overflow ----------------------------
            if P_ic_new > 1e6
                E_bwd = E_bwd .* sqrt(1e6 / P_ic_new);
            end
            E_ic = E_bwd;
        end

        P_ic = sum(abs(E_ic(:)).^2) * p.dx^2;
    end

    P_out       = (1 - R_OC) * P_ic;
    P_out_vec(m)= max(P_out, 0);
    P_ic_vec(m) = P_ic;
    f_TL_vec(m) = f_TL;
    T_max_vec(m)= max(T_3d(:));

    if verbose
        fprintf('P_out=%.3f W, P_ic=%.2f W, f_TL=%.0f mm, T_max=%.1f K  [%.1fs]\n',...
                P_out_vec(m), P_ic, f_TL*1e3, T_max_vec(m), toc);
    end
end

% ---- Threshold and slope efficiency ------------------------------------
above = P_out_vec > 1e-6;
P_thr = NaN; slope_eff = NaN;
if sum(above) >= 2
    % Linear fit above threshold
    idx_above = find(above);
    coeffs    = polyfit(P_pump_vec(idx_above), P_out_vec(idx_above), 1);
    slope_eff = coeffs(1);
    P_thr     = -coeffs(2) / coeffs(1);   % x-intercept
    if verbose
        fprintf('\nThreshold ≈ %.2f W,  slope efficiency ≈ %.1f%%\n', ...
                P_thr, slope_eff*100);
    end
end

results.P_pump         = P_pump_vec;
results.P_out          = P_out_vec;
results.P_intracavity  = P_ic_vec;
results.P_threshold    = P_thr;
results.slope_eff      = slope_eff;
results.f_thermal      = f_TL_vec;
results.T_max          = T_max_vec;
results.p              = p;

plot_cavity_results(results);

end

% =========================================================================
function E = make_gaussian_field(p, P, w0)
A0 = sqrt(2*P/(pi*w0^2));
E  = A0 .* exp(-p.R_grid.^2 ./ w0^2);
end

function plot_cavity_results(res)
p = res.p;
figure('Name',sprintf('CW Flat/Flat Laser [%s]',p.label),...
       'NumberTitle','off','Position',[140 140 1000 650]);

subplot(2,2,1);
plot(res.P_pump, res.P_out,'b-o','LineWidth',2,'MarkerSize',7);
hold on;
if ~isnan(res.P_threshold)
    xline(res.P_threshold,'r--','LineWidth',1.5,'Label',sprintf('P_{th}=%.1f W',res.P_threshold));
    % Slope line
    Pp_fit = linspace(res.P_threshold, max(res.P_pump), 50);
    plot(Pp_fit, res.slope_eff*(Pp_fit - res.P_threshold), 'r-','LineWidth',1.5);
end
xlabel('Pump power (W)'); ylabel('Output power (W)');
title('CW output power'); grid on; legend('Simulation','Threshold','Slope fit');

subplot(2,2,2);
plot(res.P_pump, res.P_intracavity,'m-s','LineWidth',2,'MarkerSize',7);
xlabel('Pump power (W)'); ylabel('Intracavity power (W)');
title('Intracavity power'); grid on;

subplot(2,2,3);
valid = isfinite(res.f_thermal) & res.f_thermal ~= 0;
if any(valid)
    plot(res.P_pump(valid), res.f_thermal(valid)*1e3,'g-d','LineWidth',2,'MarkerSize',7);
    xlabel('Pump power (W)'); ylabel('f_{thermal} (mm)');
    title('Thermal lens focal length'); grid on;
end

subplot(2,2,4);
plot(res.P_pump, res.T_max,'k-^','LineWidth',2,'MarkerSize',7);
xlabel('Pump power (W)'); ylabel('T_{max} (K)');
title('Peak crystal temperature'); grid on;

sgtitle(sprintf('CW Flat/Flat Cavity: %s', p.label),'FontWeight','bold');
drawnow;
end
