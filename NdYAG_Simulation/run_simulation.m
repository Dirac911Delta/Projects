% run_simulation.m
%
% Top-level entry point for the Nd:YAG laser simulation suite.
%
% This script runs all three scenarios based on the model of:
%   Harrison, Forbes & Naidoo, Opt. Express 31(15), 24516 (2023)
%   "Improving performance prediction of diode end-pumped solid-state
%    Nd:YAG rod amplifiers by incorporating pump mode evolution"
%
% SCENARIO 1 — Single-pump amplifier (paper baseline)
%   Reproduces the paper's experimental results using the exact crystal and
%   beam parameters from Table 1 (L=25 mm, R=2 mm, contra-propagating).
%
% SCENARIO 2 — Dual-pump amplifier (user extension)
%   Applies the same physics model to the user's larger crystal
%   (L=70 mm, R=3 mm, 0.5% Nd:YAG) with pump beams entering both end-faces.
%
% SCENARIO 3 — CW flat/flat laser oscillator (user extension)
%   Simulates a CW laser using the same crystal as Scenario 2, placed in a
%   flat/flat resonator cavity (HR + OC), using Fox-Li round-trip iteration.
%
% USAGE
%   Simply run this script in MATLAB 2023a (or later):
%       >> run_simulation
%
%   To run individual scenarios:
%       >> res1 = single_pump_amplifier();          % paper baseline
%       >> res2 = dual_pump_amplifier();            % dual pump
%       >> res3 = CW_flat_flat_laser();             % CW cavity
%
%   To use the user's crystal instead of the paper's crystal for Scenario 1:
%       >> res1u = single_pump_amplifier([], [], 'contra', 'user');
%
% All figures are generated automatically. Results are returned as structs.
% --------------------------------------------------------------------------
% Add the simulation folder to the path if not already there
script_dir = fileparts(mfilename('fullpath'));
if ~isempty(script_dir) && ~any(strcmp(path, script_dir))
    addpath(script_dir);
end

clear; close all; clc;

fprintf('============================================================\n');
fprintf('  Nd:YAG End-Pumped Laser Simulation Suite\n');
fprintf('  Based on: Harrison, Forbes & Naidoo,\n');
fprintf('  Opt. Express 31(15), 24516 (2023)\n');
fprintf('============================================================\n\n');

%% -----------------------------------------------------------------------
%  SCENARIO 1: Single-pump amplifier  (reproduces paper results)
%  Crystal: L=25 mm, R=2 mm, 0.5%% Nd:YAG
%  Pump:    MMFC diode 808 nm, DOE-shaped, contra-propagating
%  Seed:    Gaussian TEM00, 100 mW, 1064 nm
% -----------------------------------------------------------------------
fprintf('>>> Scenario 1: Single-pump amplifier (paper baseline)\n');

P_pump_paper = [0, 5, 10, 15, 20, 25, 30, 35, 38];   % Table 1 range [W]

res1 = single_pump_amplifier(P_pump_paper, 0.1, 'contra', 'paper', true);

fprintf('\n--- Scenario 1 Summary ---\n');
fprintf('P_pump (W)  |  P_out (W)  |  Gain (dB)  |  f_TL (mm)  |  T_max (K)\n');
fprintf('-----------------------------------------------------------------------\n');
for k = 1:numel(P_pump_paper)
    fprintf('%10.1f  | %11.3f | %11.1f | %11.1f | %10.1f\n', ...
            P_pump_paper(k), res1.P_out(k), res1.gain_dB(k), ...
            res1.f_thermal(k)*1e3, res1.T_max(k));
end

%% -----------------------------------------------------------------------
%  SCENARIO 2: Dual-pump amplifier  (user's extended crystal)
%  Crystal: L=70 mm, R=3 mm, 0.5%% Nd:YAG
%  Pump:    MMFC diode 808 nm, DOE-shaped, entering both end-faces
%  Seed:    Gaussian TEM00, 100 mW, 1064 nm
% -----------------------------------------------------------------------
fprintf('\n\n>>> Scenario 2: Dual-pump amplifier (user crystal L=70 mm)\n');

P_pump_dual = [0, 5, 10, 20, 30, 40, 50, 75, 100];   % total pump power [W]

res2 = dual_pump_amplifier(P_pump_dual, 0.1, 'user', true);

fprintf('\n--- Scenario 2 Summary ---\n');
fprintf('P_total (W) |  P_out (W)  |  Gain (dB)  |  f_TL (mm)  |  T_max (K)\n');
fprintf('-----------------------------------------------------------------------\n');
for k = 1:numel(P_pump_dual)
    fprintf('%10.1f  | %11.3f | %11.1f | %11.1f | %10.1f\n', ...
            P_pump_dual(k), res2.P_out(k), res2.gain_dB(k), ...
            res2.f_thermal(k)*1e3, res2.T_max(k));
end

%% -----------------------------------------------------------------------
%  SCENARIO 3: CW flat/flat laser oscillator  (user's extended crystal)
%  Crystal:  L=70 mm, R=3 mm, 0.5%% Nd:YAG
%  Cavity:   flat HR (R=0.9998) + flat OC (R=0.90), flat/flat geometry
%  Pump:     same as Scenario 2 (single-end, contra-propagating)
% -----------------------------------------------------------------------
fprintf('\n\n>>> Scenario 3: CW flat/flat laser oscillator\n');

P_pump_cw = [0, 2, 5, 8, 10, 15, 20, 25, 30, 35, 38];

res3 = CW_flat_flat_laser(P_pump_cw, 'user', true);

fprintf('\n--- Scenario 3 Summary ---\n');
if ~isnan(res3.P_threshold)
    fprintf('Threshold pump power  : %.2f W\n', res3.P_threshold);
    fprintf('Slope efficiency      : %.1f%%\n', res3.slope_eff*100);
end
fprintf('P_pump (W)  |  P_out (W)  |  P_ic (W)  |  f_TL (mm)\n');
fprintf('----------------------------------------------------------\n');
for k = 1:numel(P_pump_cw)
    fprintf('%10.1f  | %11.3f | %10.3f | %11.1f\n', ...
            P_pump_cw(k), res3.P_out(k), res3.P_intracavity(k), ...
            res3.f_thermal(k)*1e3);
end

%% -----------------------------------------------------------------------
%  Summary comparison figure
% -----------------------------------------------------------------------
figure('Name','Summary Comparison','NumberTitle','off','Position',[160 50 1100 500]);

subplot(1,3,1);
plot(res1.P_pump, res1.P_out,'b-o','LineWidth',2,'MarkerSize',7);
xlabel('Pump power (W)'); ylabel('Output power (W)');
title('Scenario 1: Single pump (paper)'); grid on;

subplot(1,3,2);
plot(res2.P_pump, res2.P_out,'r-s','LineWidth',2,'MarkerSize',7);
xlabel('Total pump power (W)'); ylabel('Output power (W)');
title('Scenario 2: Dual pump (user crystal)'); grid on;

subplot(1,3,3);
plot(res3.P_pump, res3.P_out,'m-d','LineWidth',2,'MarkerSize',7);
hold on;
if ~isnan(res3.P_threshold)
    xline(res3.P_threshold,'k--','LineWidth',1.5);
end
xlabel('Pump power (W)'); ylabel('Output power (W)');
title('Scenario 3: CW flat/flat cavity'); grid on;

sgtitle('Nd:YAG Laser Simulation — All Scenarios','FontWeight','bold','FontSize',13);
drawnow;

fprintf('\n============================================================\n');
fprintf('  Simulation complete. All figures generated.\n');
fprintf('============================================================\n');
