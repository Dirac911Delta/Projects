function material = materials(materialName)
% materials.m
% Material database for Nd-doped laser glasses.
% Data extracted from three manufacturer tables:
%   Table 1 - Optical Specifications
%   Table 2 - Laser Specifications
%   Table 3 - Thermal and Other Specifications
%
% SI conversion summary (applied below):
%   dn/dT          [1e-6 /°C]      -> multiply by 1e-6  -> [1/K]
%   emission_cross_section [1e-20 cm²] -> multiply by 1e-24 -> [m²]
%   effective_bandwidth [nm]        -> multiply by 1e-9  -> [m]
%   fluorescence_wavelength [nm]    -> multiply by 1e-9  -> [m]
%   Nd_concentration [1e20 ions/cm³]-> multiply by 1e26  -> [ions/m³]
%   thermal_conductivity [W/(m·K)]  -> already SI (table notation W/Mk = W/(m·K))
%   specific_heat [J/(g·K)]         -> multiply by 1e3   -> [J/(kg·K)]
%   density [g/cm³]                 -> multiply by 1e3   -> [kg/m³]
%   thermal_expansion [1e-7 /K]     -> multiply by 1e-7  -> [1/K]
%   thermal_opl_coeff [1e-6 /K]     -> multiply by 1e-6  -> [1/K]
%   youngs_modulus [GPa]            -> multiply by 1e9   -> [Pa]
%   lifetime [µs]                   -> multiply by 1e-6  -> [s]
%   absorption_coeff [cm⁻¹]         -> multiply by 1e2   -> [m⁻¹]
%
% Missing table values (\) are stored as NaN.
% absorption_808 is NOT in any table; it remains NaN and must be
% user-supplied or estimated (see pump_model.m / main.m).

db = build_database();
idx = find(strcmpi({db.name}, materialName), 1);
if isempty(idx)
    error('Unknown material "%s". Valid names: %s', materialName, strjoin({db.name}, ' '));
end
material = db(idx);
end

% =========================================================================
function db = build_database()

% --- Table 1: Optical Specifications ---
% Nonlinear refractive index n² upper bound (x1e-13 e.s.u.)
%              N31   N41   N51   NAP2  NAP4  NF1   NF2   NSG2
n2_esu =      [1.2,  1.04, 1.04, 1.25, 1.10, 0.6,  0.86, 1.6 ];

% Refractive index at 1053 nm (central value ± 0.003 or ± 0.005)
n1053 =       [1.535,1.504,1.505,1.537,1.515,1.464,1.514,1.560];

% Abbe value (dimensionless)
abbe =        [65.6, 68.2, 68.2, 67,   67,   88,   77,   59  ];

% dn/dT (1e-6 /°C, range 20–100°C) -> stored as [1/K]
% N41 is missing (\)
dndT_raw =    [-4.3, NaN,  -9.0, -9.0, 1.9, -8.84,-8.6,  2.0 ];
dndT =        dndT_raw * 1e-6;      % [1/K]

% --- Table 2: Laser Specifications ---
% Nominal wt% dopant (reference doping level, not user-set)
%   NaN where table shows "\"
wt_pct_nominal = [3.5, 4.6, 4.0, NaN, NaN, 0.88, 1.07, NaN];   % see notes below

% Nd³⁺ concentration (1e20 ions/cm³) -> [ions/m³]
% NAP2, NAP4, NSG2 are missing
nd_conc_raw = [3.4, 4.3, 3.9, NaN, NaN, 0.2, 1.2, NaN];   % [1e20 ions/cm³]
nd_conc =     nd_conc_raw * 1e26;                            % [ions/m³]

% Stimulated emission cross section (1e-20 cm²) -> [m²]
sigma_raw =   [3.8, 3.9, 4.3, 3.6, 3.1, 2.7, 3.4, 2.7];   % [1e-20 cm²]
sigma =       sigma_raw * 1e-24;                             % [m²]

% Effective bandwidth (nm) -> [m]
bw_nm =       [25.4, 25.5, 24.5, 25.4, 28.5, 32.8, 30.4, 34];
bw =          bw_nm * 1e-9;                                  % [m]

% Fluorescence peak wavelength (nm) -> [m]
lfl_nm =      [1053, 1053, 1053, 1052, 1052, 1053, 1053, 1060];
lfl =         lfl_nm * 1e-9;                                 % [m]

% Absorption coefficients (upper bounds, cm⁻¹) -> [m⁻¹]
% Three wavelength entries per material: 1053 nm, 400 nm, 3333 nm.
% NSG2 only has 1053 nm value listed.
%                         N31      N41      N51      NAP2     NAP4     NF1      NF2      NSG2
abs_1053_cm =           [0.0015,  0.0015,  0.0015,  0.0015,  0.002,   0.001,   0.001,   0.0015];
abs_400_cm =            [0.25,    0.25,    0.25,    0.25,    0.3,     0.04,    0.04,    NaN   ];
abs_3333_cm =           [1.5,     1.5,     1.5,     1.5,     1.5,     0.08,    0.08,    NaN   ];
abs_1053 = abs_1053_cm * 1e2;    % [m⁻¹]
abs_400  = abs_400_cm  * 1e2;    % [m⁻¹]
abs_3333 = abs_3333_cm * 1e2;    % [m⁻¹]

% NOTE: 808 nm absorption is NOT present in any table.
% It must be supplied by the user or estimated (see pump_model.m).

% --- Table 3: Thermal Specifications ---
%                         N31   N41   N51   NAP2  NAP4  NF1   NF2   NSG2
T_transition_C =         [445,  467,  408,  500,  545,  450,  490,  485 ];  % [°C]
T_softening_C =          [485,  503,  448,  550,  600,  491,  528,  530 ];  % [°C]

% Linear thermal expansion coefficient (1e-7 /K, 30–100°C) -> [1/K]
cte_raw =                [116,  129,  141,  87,   63,   152,  142,  95  ];  % [1e-7 /K]
cte =                    cte_raw * 1e-7;                         % [1/K]

% Thermal coefficient of optical path length (1e-6 /K, 50–100°C) -> [1/K]
% N41 is missing
thermal_opl_raw =        [1.4,  NaN, -1.9,  3.8,  5.0, -1.86,-1.2,  7.0 ];  % [1e-6 /K]
thermal_opl =            thermal_opl_raw * 1e-6;                 % [1/K]

% Thermal conductivity at 25°C [W/(m·K)] (table notation W/Mk)
% N41, N51, NF2 are missing
kth =                    [0.59, NaN,  NaN,  0.76, 0.88, 0.865,NaN,  1.2 ];  % [W/(m·K)]

% Specific heat capacity at 25°C [J/(g·K)] -> [J/(kg·K)]
% N41, N51, NF1, NF2, NSG2 are missing
cp_raw =                 [0.75, NaN,  NaN,  0.757,0.775,NaN,  NaN,  NaN ];  % [J/(g·K)]
cp =                     cp_raw * 1e3;                           % [J/(kg·K)]

% --- Table 3: Other Specifications ---
% Density [g/cm³] -> [kg/m³]
% NSG2 is not listed in this table
rho_raw =                [2.87, 2.62, 2.7,  2.84, 2.58, 3.65, 3.68, NaN];  % [g/cm³]
rho =                    rho_raw * 1e3;                          % [kg/m³]

% Young's modulus [GPa] -> [Pa]
% NSG2 is not listed
E_raw =                  [58.3, 52.4, 45.2, 58,   67,   73,   76,   NaN];  % [GPa]
E =                      E_raw * 1e9;                            % [Pa]

% Poisson's ratio (dimensionless)
% NF1, NF2, NSG2 are missing
nu =                     [0.26, 0.25, 0.26, 0.25, 0.25, NaN,  NaN,  NaN];

% Knoop hardness [kg/cm²]
% NSG2 not listed
knoop =                  [404,  347,  302,  382,  549,  343,  423,  NaN];

% Fracture toughness [MPa·m^(1/2)]
% NSG2 not listed
KIC =                    [0.58, 0.62, 0.66, 0.68, 0.74, 0.35, 0.58, NaN];  % [MPa·m^0.5]

% Deliquescence coefficient (H₂O 98°C) [mg/(cm²/day)]
% N31, NF1, NF2, NSG2 are missing
deliq =                  [NaN,  0.41, 2.2,  0.003,0.002,NaN,  NaN,  NaN];

% =========================================================================
% Lifetime-vs-doping tables (µs -> s at specified wt%)
% Values are guaranteed minima (≥); stored as the minimum guaranteed value.
% N31, N41, N51, NAP2, NAP4, NSG2 use Nd₂O₃ wt%.
% NF1, NF2 use NdF₃ wt%.
% =========================================================================
lifetime_tables = {
% N31  doping [Nd2O3 wt%]        lifetime [µs]
    [0.5, 1.2, 3.5, 4.2],        [370, 360, 315, 310];
% N41
    [0.5, 1.2, 3.5, 4.6],        [370, 360, 315, 310];
% N51
    [0.5, 1.2, 3.5, 4.2],        [375, 365, 320, 315];
% NAP2
    [0.5, 1.0, 2.0, 3.0],        [360, 350, 330, 310];
% NAP4
    [0.5, 1.0, 2.0, 3.0],        [370, 360, 330, 310];
% NF1  doping [NdF3 wt%]
    [0.53, 1.07],                 [515, 495];
% NF2
    [0.53, 1.07],                 [430, 410];
% NSG2  doping [Nd2O3 wt%]
    [0.5, 1.0, 2.0, 3.0],        [380, 360, 330, 270];
};

names = {'N31','N41','N51','NAP2','NAP4','NF1','NF2','NSG2'};
n = numel(names);

% Allocate struct array
db = repmat(struct( ...
    'name',                        '', ...
    'refractive_index_1053nm',     NaN, ...
    'dn_dT',                       NaN, ...
    'nonlinear_n2_esu',            NaN, ...
    'abbe_value',                  NaN, ...
    'emission_cross_section',      NaN, ...
    'effective_bandwidth',         NaN, ...
    'fluorescence_wavelength',     NaN, ...
    'wt_pct_nominal',              NaN, ...
    'Nd_concentration',            NaN, ...
    'lifetime_vs_doping',          struct(), ...
    'absorption_808',              NaN, ...
    'absorption_coeff_1053nm_max', NaN, ...
    'absorption_coeff_400nm_max',  NaN, ...
    'absorption_coeff_3333nm_max', NaN, ...
    'thermal_conductivity',        NaN, ...
    'specific_heat',               NaN, ...
    'density',                     NaN, ...
    'thermal_expansion_coeff',     NaN, ...
    'thermal_opl_coeff',           NaN, ...
    'transition_temperature_C',    NaN, ...
    'softening_temperature_C',     NaN, ...
    'youngs_modulus_Pa',           NaN, ...
    'poissons_ratio',              NaN, ...
    'knoop_hardness_kg_per_cm2',   NaN, ...
    'fracture_toughness_MPa_sqrtm',NaN, ...
    'deliquescence_mg_per_cm2_day',NaN), 1, n);

for k = 1:n
    db(k).name                         = names{k};
    db(k).refractive_index_1053nm      = n1053(k);
    db(k).dn_dT                        = dndT(k);          % [1/K]
    db(k).nonlinear_n2_esu             = n2_esu(k);        % [1e-13 e.s.u.], upper bound
    db(k).abbe_value                   = abbe(k);
    db(k).emission_cross_section       = sigma(k);         % [m²]
    db(k).effective_bandwidth          = bw(k);            % [m]
    db(k).fluorescence_wavelength      = lfl(k);           % [m]
    db(k).wt_pct_nominal               = wt_pct_nominal(k);% [wt%] nominal dopant compound
    db(k).Nd_concentration             = nd_conc(k);       % [ions/m³]
    db(k).absorption_808               = NaN;              % NOT IN TABLES - user must supply
    db(k).absorption_coeff_1053nm_max  = abs_1053(k);      % [m⁻¹] upper bound
    db(k).absorption_coeff_400nm_max   = abs_400(k);       % [m⁻¹] upper bound
    db(k).absorption_coeff_3333nm_max  = abs_3333(k);      % [m⁻¹] upper bound
    db(k).thermal_conductivity         = kth(k);           % [W/(m·K)]
    db(k).specific_heat                = cp(k);            % [J/(kg·K)]
    db(k).density                      = rho(k);           % [kg/m³]
    db(k).thermal_expansion_coeff      = cte(k);           % [1/K]
    db(k).thermal_opl_coeff            = thermal_opl(k);   % [1/K]
    db(k).transition_temperature_C     = T_transition_C(k);% [°C]
    db(k).softening_temperature_C      = T_softening_C(k); % [°C]
    db(k).youngs_modulus_Pa            = E(k);             % [Pa]
    db(k).poissons_ratio               = nu(k);
    db(k).knoop_hardness_kg_per_cm2    = knoop(k);
    db(k).fracture_toughness_MPa_sqrtm = KIC(k);
    db(k).deliquescence_mg_per_cm2_day = deliq(k);

    % Lifetime lookup table with interpolation closure
    doping_table  = lifetime_tables{k, 1};
    lifetime_table = lifetime_tables{k, 2} * 1e-6;   % µs -> s
    db(k).lifetime_vs_doping = struct( ...
        'doping_wt_pct', doping_table, ...
        'lifetime_s',    lifetime_table, ...
        'interp',        @(x)interp_lifetime(doping_table, lifetime_table, x));
end
end

% =========================================================================
function tau = interp_lifetime(doping_wt_pct, lifetime_s, query)
% interp_lifetime  Interpolate fluorescence lifetime vs Nd dopant wt%.
% Values outside the measured range return NaN with a warning.
% Inputs are guaranteed minima from table; interpolant passes through them.
valid = ~(isnan(doping_wt_pct) | isnan(lifetime_s));
if sum(valid) == 0
    tau = NaN(size(query));
    return;
elseif sum(valid) == 1
    tau = lifetime_s(valid) * ones(size(query));
    return;
end
x = doping_wt_pct(valid);
y = lifetime_s(valid);
tau = interp1(x, y, query, 'pchip', NaN);
out_values = query(query < min(x) | query > max(x));
if ~isempty(out_values)
    warning('Lifetime query outside measured range [%.4g, %.4g] wt%%. Out-of-range values return NaN. Queried: %s', ...
        min(x), max(x), mat2str(out_values));
end
end
