function material = materials(materialName)
% materials.m
% Material database for Nd-doped laser glasses from provided tables.
% Missing table values are stored as NaN.
% SI conversion notes:
%   dn/dT values in source are [1e-6 / degC] -> multiplied by 1e-6 [1/K].

db = build_database();
idx = find(strcmpi({db.name}, materialName), 1);
if isempty(idx)
    error('Unknown material "%s".', materialName);
end
material = db(idx);
end

function db = build_database()
names = {'N31','N41','N51','NAP2','NAP4','NF1','NF2','NSG2'};

% Extracted from provided optical table image (where visible)
n1053 = [1.535, 1.504, 1.505, 1.537, 1.515, 1.464, 1.514, 1.560];
dndT_micro_per_K = [-4.3, NaN, -9.0, -9.0, 1.9, -8.84, -8.6, 2.0];

% Values not visible in available image set are NaN by rule.
nanv = NaN(1, numel(names));

db = repmat(struct( ...
    'name', '', ...
    'refractive_index_1053nm', NaN, ...
    'dn_dT', NaN, ...
    'emission_cross_section', NaN, ...
    'effective_bandwidth', NaN, ...
    'fluorescence_wavelength', NaN, ...
    'Nd_concentration', NaN, ...
    'lifetime_vs_doping', struct('doping_wt_pct', NaN, 'lifetime_s', NaN, 'interp', @(x) nan(size(x))), ...
    'absorption_808', NaN, ...
    'thermal_conductivity', NaN, ...
    'specific_heat', NaN, ...
    'density', NaN), 1, numel(names));

for k = 1:numel(names)
    db(k).name = names{k};
    db(k).refractive_index_1053nm = n1053(k);
    db(k).dn_dT = dndT_micro_per_K(k) * 1e-6;  % [1/K]
    db(k).emission_cross_section = nanv(k);      % [m^2]
    db(k).effective_bandwidth = nanv(k);         % [m] or [Hz], user-defined convention
    db(k).fluorescence_wavelength = nanv(k);     % [m]
    db(k).Nd_concentration = nanv(k);            % user-provided basis
    db(k).absorption_808 = nanv(k);              % [1/m]

    % Lifetime-vs-doping placeholder; user should populate from tabulated data.
    doping = NaN;
    lifetime = NaN;
    db(k).lifetime_vs_doping = struct( ...
        'doping_wt_pct', doping, ...
        'lifetime_s', lifetime, ...
        'interp', @(x)interp_lifetime(doping, lifetime, x));
end
end

function tau = interp_lifetime(doping_wt_pct, lifetime_s, query)
% Robust interpolation helper allowing missing/partial datasets.
if all(isnan(doping_wt_pct)) || all(isnan(lifetime_s))
    tau = NaN(size(query));
    return;
end
valid = ~(isnan(doping_wt_pct) | isnan(lifetime_s));
if sum(valid) == 0
    tau = NaN(size(query));
elseif sum(valid) == 1
    tau = lifetime_s(valid) * ones(size(query));
else
    tau = interp1(doping_wt_pct(valid), lifetime_s(valid), query, 'pchip', 'extrap');
end
end
