function results = heat_solver(material, pump, settings)
% heat_solver.m
% Axisymmetric transient heat equation in (r,z):
%   dT/dt = a * [ (1/r)*d/dr(r dT/dr) + d2T/dz2 ] + Q/(rho*cp)
% Crank-Nicolson time integration:
%   (I - 0.5*dt*a*L) T^{n+1} = (I + 0.5*dt*a*L) T^n + dt*S

r = pump.r;
z = pump.z;
Q = pump.Q;

nr = numel(r);
nz = numel(z);
N = nr * nz;
dr = r(2) - r(1);
dz = z(2) - z(1);

dt = settings.dt_s;
t = 0:dt:settings.t_end_s;
nt = numel(t);

k = material.thermal_conductivity;
rho = material.density;
cp = material.specific_heat;
a = k / (rho * cp);

S = Q / (rho * cp);  % [K/s]

% Sparse Laplacian operator L
Lop = spalloc(N, N, 7*N);

for j = 1:nz
    for i = 1:nr
        p = idx(i, j, nr);

        at_r_outer = (i == nr);
        at_z0 = (j == 1);
        at_zL = (j == nz);

        % Dirichlet rows handled later in A/B; keep Laplacian row zero
        if at_r_outer
            continue;
        end
        if at_z0 && strcmpi(settings.bc_z0.type, 'fixed')
            continue;
        end
        if at_zL && strcmpi(settings.bc_zL.type, 'fixed')
            continue;
        end

        if i == 1
            % r=0 symmetry: radial operator limit
            Lop(p, p) = Lop(p, p) - 4/dr^2;
            Lop(p, idx(i+1,j,nr)) = Lop(p, idx(i+1,j,nr)) + 4/dr^2;
        else
            ri = r(i);
            ar = 1/dr^2;
            br = 1/(2*ri*dr);
            Lop(p, idx(i-1,j,nr)) = Lop(p, idx(i-1,j,nr)) + (ar - br);
            Lop(p, p)            = Lop(p, p) - 2*ar;
            Lop(p, idx(i+1,j,nr)) = Lop(p, idx(i+1,j,nr)) + (ar + br);
        end

        if j == 1
            % z=0 insulated
            Lop(p, p) = Lop(p, p) - 2/dz^2;
            Lop(p, idx(i,j+1,nr)) = Lop(p, idx(i,j+1,nr)) + 2/dz^2;
        elseif j == nz
            % z=L insulated
            Lop(p, p) = Lop(p, p) - 2/dz^2;
            Lop(p, idx(i,j-1,nr)) = Lop(p, idx(i,j-1,nr)) + 2/dz^2;
        else
            Lop(p, idx(i,j-1,nr)) = Lop(p, idx(i,j-1,nr)) + 1/dz^2;
            Lop(p, p)            = Lop(p, p) - 2/dz^2;
            Lop(p, idx(i,j+1,nr)) = Lop(p, idx(i,j+1,nr)) + 1/dz^2;
        end
    end
end

I = speye(N);
A = I - 0.5 * dt * a * Lop;
B = I + 0.5 * dt * a * Lop;

% Initial condition: heat sink temperature
T = settings.heat_sink_temperature_K * ones(N, 1);
Svec = reshape(S, [], 1);

Tall = zeros(nr, nz, nt);
Tall(:,:,1) = reshape(T, nr, nz);

% Precompute fixed-temperature row masks and values
[fixedMask, fixedValue] = build_fixed_mask(settings, nr, nz);
fixedIdx = find(fixedMask);

for n = 2:nt
    rhs = B * T + dt * Svec;

    % Enforce fixed boundaries in linear system
    Aeff = A;
    if ~isempty(fixedIdx)
        Aeff(fixedIdx, :) = 0;
        Aeff(fixedIdx, fixedIdx) = speye(numel(fixedIdx));
        rhs(fixedIdx) = fixedValue(fixedIdx);
    end

    T = Aeff \ rhs;
    Tall(:,:,n) = reshape(T, nr, nz);
end

results = struct();
results.r = r;
results.z = z;
results.t = t;
results.T = Tall;
results.Q = Q;
end

function p = idx(i, j, nr)
p = i + (j - 1) * nr;
end

function [mask, val] = build_fixed_mask(settings, nr, nz)
N = nr * nz;
mask = false(N,1);
val = settings.heat_sink_temperature_K * ones(N,1);

% r = R fixed heat sink
for j = 1:nz
    p = nr + (j-1)*nr;
    mask(p) = true;
    val(p) = settings.heat_sink_temperature_K;
end

% z boundaries configurable
if strcmpi(settings.bc_z0.type, 'fixed')
    for i = 1:nr
        p = i;
        mask(p) = true;
        val(p) = settings.bc_z0.value_K;
    end
end

if strcmpi(settings.bc_zL.type, 'fixed')
    for i = 1:nr
        p = i + (nz-1)*nr;
        mask(p) = true;
        val(p) = settings.bc_zL.value_K;
    end
end
end
