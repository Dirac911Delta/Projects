function [E_out, I_abs] = bpm_propagate(E_in, p, dz, lambda0, n_med, ...
                                         alpha_2d, delta_n_2d, gain_2d)
% BPM_PROPAGATE  One axial step of the split-step angular-spectrum BPM.
%
% Implements the symmetric split-step scheme described in Harrison et al.,
% Opt. Express 31(15), 24516 (2023), Section 3.5 (iterative Fourier beam
% propagation method).  At each step dz the algorithm:
%   1. Applies a half-step thin-lens mask  (absorption + gain + phase)
%   2. Propagates the field in free space by dz via the angular-spectrum
%      transfer function  H = exp(i·kz·dz)
%   3. Applies the second half-step thin-lens mask
%
% The angular-spectrum transfer function uses MATLAB's un-shifted FFT
% convention so no fftshift is needed on the field arrays.
%
% INPUTS
%   E_in      [N×N] complex field amplitude [sqrt(W)/m],  |E|² = I [W/m²]
%   p               parameter struct (NdYAG_params)
%   dz              propagation step                        [m]
%   lambda0         free-space wavelength                   [m]
%   n_med           background refractive index (scalar)
%   alpha_2d  [N×N] power absorption coefficient           [m⁻¹]  (or scalar)
%   delta_n_2d[N×N] refractive-index perturbation Δn       (or scalar, def=0)
%   gain_2d   [N×N] power gain coefficient g               [m⁻¹]  (or scalar, def=0)
%
% OUTPUTS
%   E_out     [N×N] complex field after one step
%   I_abs     [N×N] absorbed intensity in this slice        [W/m²]
%
% NOTE: amplitude absorption  → field multiplied by exp(-α/2 · dz)
%       amplitude gain        → field multiplied by exp(+g/2 · dz)
%       refractive-index      → field multiplied by exp(+i·k₀·Δn · dz)

if nargin < 7 || isempty(delta_n_2d), delta_n_2d = 0; end
if nargin < 8 || isempty(gain_2d),    gain_2d    = 0; end

N  = p.N_grid;
dx = p.dx;

k0 = 2*pi / lambda0;          % free-space wavenumber  [m⁻¹]
k  = k0 * n_med;               % wavenumber in medium

% ---- Angular-spectrum transfer function ---------------------------------
% Spatial-frequency grid (un-shifted FFT convention: DC at index 1)
dkx    = 2*pi / (N * dx);
kx_vec = [0 : N/2-1, -N/2 : -1] * dkx;    % row vector [1×N]
[KX, KY] = meshgrid(kx_vec, kx_vec);       % [N×N]

KZ2 = k^2 - KX.^2 - KY.^2;
KZ  = sqrt(max(KZ2, 0));       % real kz only (evanescent terms zeroed)
H   = exp(1i .* KZ .* dz);    % free-space transfer function [N×N]

% ---- Half-step thin-lens mask -------------------------------------------
%   amplitude factor: exp((-α+g) · |dz|/2) applied twice = exp((-α+g) · |dz|)
%     Uses abs(dz) so gain and absorption are always positive-definite
%     regardless of propagation direction (dz may be negative for backward pass).
%   phase factor: exp(i·k₀·Δn · dz/2) applied twice = exp(i·k₀·Δn · dz)
%     Uses signed dz so that accumulated phase is correctly reversed for
%     backward propagation.
amp_half   = exp((-alpha_2d + gain_2d) .* (abs(dz)/2));  % [N×N] or scalar
phase_half = exp(1i .* k0 .* delta_n_2d .* (dz/2));      % [N×N] or scalar
mask_half  = amp_half .* phase_half;

% ---- Step 1: first half-step thin-lens mask -----------------------------
E1 = E_in .* mask_half;

% ---- Step 2: free-space propagation (angular-spectrum method) -----------
E_ft   = fft2(E1);
E_prop = E_ft .* H;
E2     = ifft2(E_prop);

% ---- Step 3: second half-step thin-lens mask ----------------------------
E_out = E2 .* mask_half;

% ---- Absorbed intensity (local, in this slice) --------------------------
% dI/dz = -alpha · I  →  I_abs ≈ alpha · |E_in|² · dz  (per unit volume)
% We report total absorbed intensity over the slice (power per unit area):
I_abs = max(abs(E_in).^2 - abs(E_out).^2, 0);

end
