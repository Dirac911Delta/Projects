function p = NdYAG_params(crystal_config)
% NdYAG_params  Return a struct of physical parameters for the Nd:YAG
%               laser simulation.
%
% Usage:
%   p = NdYAG_params()              % paper crystal  (L=25 mm, R=2 mm)
%   p = NdYAG_params('paper')       % same as above
%   p = NdYAG_params('user')        % user crystal   (L=70 mm, R=3 mm)
%
% All crystal properties are taken directly from Table 1 of:
%   Harrison, Forbes & Naidoo, Opt. Express 31(15), 24516 (2023)
%   "Improving performance prediction of diode end-pumped solid-state
%    Nd:YAG rod amplifiers by incorporating pump mode evolution"
%
% All quantities in SI units (m, W, K, s, J) unless explicitly noted.
% --------------------------------------------------------------------------

if nargin < 1, crystal_config = 'paper'; end

%% -----------------------------------------------------------------------
%  Physical constants
%  -----------------------------------------------------------------------
p.h  = 6.62607015e-34;   % Planck constant               [J·s]
p.c  = 2.99792458e8;     % Speed of light in vacuum       [m/s]
p.kB = 1.380649e-23;     % Boltzmann constant             [J/K]
p.eps0 = 8.8541878e-12;  % Permittivity of free space     [F/m]

%% -----------------------------------------------------------------------
%  Crystal geometry  (Table 1, Harrison et al. 2023)
%  -----------------------------------------------------------------------
switch lower(crystal_config)
    case 'paper'
        % ---- Paper validation crystal (Section 2, Table 1) --------------
        p.L_c   = 25e-3;    % Crystal length                [m]
        p.R_c   = 2e-3;     % Crystal radius                [m]  (4 mm diam)
        p.label = 'Paper crystal (L=25 mm, R=2 mm)';
    case 'user'
        % ---- User-specified crystal (problem statement) -----------------
        p.L_c   = 70e-3;    % Crystal length                [m]
        p.R_c   = 3e-3;     % Crystal radius                [m]  (6 mm diam)
        p.label = 'User crystal (L=70 mm, R=3 mm)';
    otherwise
        error('NdYAG_params: unknown crystal_config "%s". Use "paper" or "user".',...
              crystal_config);
end

%% -----------------------------------------------------------------------
%  Nd:YAG material properties  (Table 1, Harrison et al. 2023)
%  -----------------------------------------------------------------------
% ----- Optical / spectroscopic ------------------------------------------
p.n0            = 1.82;           % Refractive index at 1064 nm
p.doping_at     = 0.5;            % Nd doping level           [at.%]
p.N_tot         = 0.69e20 * 1e6;  % Active ion density        [m⁻³]
%                                    (0.69×10²⁰ cm⁻³ from Table 1)

p.tau_em        = 240e-6;         % Fluorescence lifetime     [s]
%                                    (τ_em = 240 μs, Table 1)

% Emission cross-section at 1064 nm (σ^s_em, Table 1)
p.sigma_em      = 2.8e-19 * 1e-4; % = 2.8e-23 m²
%                                    (2.8×10⁻¹⁹ cm², Table 1)

% Peak absorption cross-section at ~808 nm (σ^peak_abs, Table 1)
p.sigma_abs_peak = 5e-20 * 1e-4;  % = 5e-24 m²
%                                    (5×10⁻²⁰ cm², Table 1)

% Dynamic absorption cross-section range (σ^p_abs, Table 1)
% This is the spectrally-weighted effective cross-section for the
% broadband MMFC diode pump, calibrated from the measured emission spectrum.
% Range 5.9–9.1×10⁻²¹ cm² depending on diode drive current / temperature.
% The simulation uses a single weighted value; adjust via set_pump_spectrum.
p.sigma_abs_dyn_lo = 5.9e-21 * 1e-4;  % Lower bound  [m²]
p.sigma_abs_dyn_hi = 9.1e-21 * 1e-4;  % Upper bound  [m²]
p.sigma_abs_dyn    = 7.0e-21 * 1e-4;  % Default mid-range value [m²]
%                                        (≈ 7×10⁻²⁵ m²)

% Saturation intensity at 1064 nm (4-level, no GSA)
%   I_sat = hν_s / (σ_em · τ_em)   [W/m²]
p.lambda_s = 1064e-9;                          % Signal wavelength [m]
p.hnu_s    = p.h * p.c / p.lambda_s;           % Signal photon energy [J]
p.I_sat    = p.hnu_s / (p.sigma_em * p.tau_em); % ≈ 2.9×10⁷ W/m²

% ----- Mechanical / thermo-mechanical (Table 1) -------------------------
p.E_young    = 3e10;       % Young's modulus M             [Pa]
p.nu_poisson = 0.25;       % Poisson's ratio ν             [ ]
p.alpha_th   = 7.5e-6;     % Thermal expansion coefficient  [K⁻¹] (typical YAG)

% Elasto-optic (photoelastic) coefficients for [1 1 1]-cut Nd:YAG (Table 1)
p.p11 = -0.029;            % p₁₁ elasto-optic coeff.
p.p12 =  0.0091;           % p₁₂ elasto-optic coeff.
p.p44 = -0.0615;           % p₄₄ elasto-optic coeff.

% ----- Thermal (Table 1) ------------------------------------------------
p.T_sink    = 294;         % Heat-sink temperature         [K]
p.k0_therm  = 10.5;        % Thermal conductivity at T₀    [W/(m·K)]
p.T0_therm  = 300;         % Reference temperature T₀      [K]
p.xi_therm  = -0.77;       % Thermal conductivity power-law exponent ξ
%                            k(T) = k0 · (T/T0)^ξ  (Aggarwal et al.)
p.h_c       = 1.6e4;       % Thermal conductance (surface) [W/(m²·K)]
%                            (h_c = 1.6×10⁴ W/m²/K, Table 1)

% Thermo-optic coefficient dn/dT for Nd:YAG at 1064 nm
p.dndT = 7.3e-6;           % [K⁻¹]  (standard YAG value)

% Fitting parameter for T-dependent thermal conductivity from Table 1
% Used together with k0_therm, T0_therm, xi_therm above.

%% -----------------------------------------------------------------------
%  Pump beam parameters  (Section 2, Harrison et al. 2023)
%  -----------------------------------------------------------------------
p.lambda_p   = 808e-9;     % Central pump wavelength       [m]  (802–808 nm)
p.lambda_p_lo= 802e-9;     % Lower diode wavelength        [m]
p.lambda_p_hi= 808e-9;     % Upper diode wavelength        [m]
p.delta_lambda_p = 1.85e-9;% Half-width of diode spectrum  [m]  (±1.85 nm)
p.hnu_p      = p.h * p.c / p.lambda_p;  % Pump photon energy  [J]

% Pump beam geometry (Fig. 1, Section 2)
p.omega_i    = 200e-6;     % Pump Gaussian beam radius at fiber output [m]
%                            (ω₀ = 200 μm at z = 0, face A)
p.omega_f    = 1.7e-3;     % Desired FT beam radius        [m]
%                            (ω_f = 1.7 mm at z = L, face B)
p.f_lens_DOE = 100e-3;     % Fourier lens focal length f₂  [m]
%                            (f₂ = 100 mm in 4f setup)
p.P_pump_max = 38;         % Maximum pump power            [W]

% Quantum / Stokes heat fraction
p.eta_h = 1.0 - p.lambda_p / p.lambda_s;   % ≈ 0.2406

%% -----------------------------------------------------------------------
%  Seed beam parameters  (Section 2, Table 1)
%  -----------------------------------------------------------------------
p.P_seed      = 0.100;     % Seed input power P_s          [W]
p.omega_seed_B= 181e-6;    % Seed beam radius at face B    [m]  (ω_s = 181 μm)
p.omega_seed_A= 175e-6;    % Seed waist at face A (z=0)    [m]  (ω⁰_s = 175 μm)
p.f_seed_lens = 200e-3;    % Seed focusing lens f₁         [m]
p.M2_seed     = 1.04;      % Seed beam quality M²

%% -----------------------------------------------------------------------
%  CW flat/flat cavity (Scenario 3 extension)
%  -----------------------------------------------------------------------
p.R_OC           = 0.90;   % Output-coupler reflectivity
p.R_HR           = 0.9998; % High-reflector reflectivity
p.L_cav_ext      = 0.05;   % External cavity arm length (each side) [m]
p.loss_rt_other  = 0.01;   % Additional round-trip power loss

%% -----------------------------------------------------------------------
%  Simulation grid
%  -----------------------------------------------------------------------
p.N_grid  = 256;            % Transverse grid points per side (N × N)
p.L_grid  = 3 * p.R_c * 2; % Transverse grid full width = 3× crystal diam [m]
p.N_z     = 500;            % Number of axial propagation slices
p.T_ambient = 293.15;       % Ambient temperature (for boundary cond.)  [K]

% ---- Derived grid quantities -------------------------------------------
p.dx  = p.L_grid / p.N_grid;   % Transverse pixel size                 [m]
p.dz  = p.L_c    / p.N_z;      % Axial step                            [m]

% Coordinate axes: centred on (0,0)
half    = p.N_grid / 2;
idx     = (-half : half-1)';          % column index vector
p.x_vec = idx  * p.dx;               % x-axis                            [m]
p.y_vec = p.x_vec;                   % y-axis (square grid)
[p.X, p.Y] = meshgrid(p.x_vec, p.y_vec);
p.R_grid   = sqrt(p.X.^2 + p.Y.^2); % Radial coordinate map             [m]
p.z_vec    = ((0.5 : p.N_z)') * p.dz; % Axial cell-centre positions      [m]

% Crystal aperture mask (unity inside crystal, zero outside)
p.mask_crystal = double(p.R_grid <= p.R_c);

% Wavenumbers at pump and signal wavelengths (in crystal)
p.k_s = 2*pi * p.n0 / p.lambda_s;   % Signal wavenumber in medium       [m⁻¹]
p.k_p = 2*pi * p.n0 / p.lambda_p;   % Pump   wavenumber in medium       [m⁻¹]

end
