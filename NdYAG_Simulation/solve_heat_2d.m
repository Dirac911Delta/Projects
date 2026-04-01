function [T, u_kirchhoff] = solve_heat_2d(p, Q_vol_3d, alpha_eff)
% SOLVE_HEAT_2D  Steady-state 3-D temperature distribution in the Nd:YAG
%   crystal rod using a slice-by-slice 2-D finite-difference method with
%   the Kirchhoff transform to handle the temperature-dependent thermal
%   conductivity k(T) = k₀·(T/T₀)^ξ.
%
% Physics (Salinas-Alvarado et al., ref [17] in Harrison 2023):
%
%   Steady-state heat equation (cylindrical geometry, axial symmetry):
%       ∇·[k(T)·∇T] = −Q(r,z)
%
%   Kirchhoff transform (linearises the equation):
%       u(T) = ∫_{T_sink}^{T} k(T') dT'
%            = k₀·T₀/(1+ξ) · [(T/T₀)^(1+ξ) − (T_sink/T₀)^(1+ξ)]
%
%       ∇²u = −Q(r,z)   (now linear in u)
%
%   Boundary conditions:
%       At crystal surface (r = R_c):
%           k·dT/dr|_{r=R_c} = −h_c·(T − T_sink)
%       → in terms of u:   du/dr|_{r=R_c} = −h_c·(T − T_sink)
%           (Newton's cooling: h_c from Table 1, Harrison 2023)
%       At end faces (z = 0 and z = L):
%           dT/dz = 0  (thermally insulated end faces)
%
%   After solving for u on the 2-D (x,y) grid at each z-slice, T is
%   recovered via the inverse Kirchhoff transform:
%       T(r,z) = T₀·[ u(r,z)·(1+ξ)/(k₀·T₀) + (T_sink/T₀)^(1+ξ) ]^{1/(1+ξ)}
%
% INPUTS
%   p          - parameter struct (NdYAG_params)
%   Q_vol_3d   - [N×N×Nz] volumetric heat source  Q(x,y,z)  [W/m³]
%                computed by the caller as  Q = eta_h · alpha_eff · I_p
%   alpha_eff  - effective pump absorption coeff.             [m⁻¹]  (scalar)
%                (used only for reporting; not needed inside this function)
%
% OUTPUTS
%   T          - [N×N×Nz] temperature distribution            [K]
%   u_kirchhoff- [N×N×Nz] Kirchhoff-transformed temperature  [W/m]
%
% NOTE: The 2-D FDM is solved on the (x,y) grid (N×N) at each of the Nz
%       axial slices independently (axial conduction neglected – valid when
%       the pump absorption depth >> crystal radius, Lc/R >> 1).  For the
%       25 mm × 2 mm crystal (L/R = 12.5) this is a good approximation.

N     = p.N_grid;
Nz    = p.N_z;
dx    = p.dx;
k0    = p.k0_therm;
T0    = p.T0_therm;
xi    = p.xi_therm;
T_s   = p.T_sink;
h_c   = p.h_c;
R_c   = p.R_c;
R_g   = p.R_grid;          % [N×N] radial positions

% ---- Kirchhoff reference value at the sink temperature -----------------
% The Kirchhoff transform integrates k(T) from T_sink to T:
%   u(T) = ∫_{T_sink}^{T} k(T') dT'
%        = k0·T0/(1+ξ) · [ (T/T0)^(1+ξ) − (T_sink/T0)^(1+ξ) ]
% so u = 0 at T = T_sink.  The FDM solves ∇²u = −Q with u = 0 at r = R_c.
coeff = k0 * T0 / (1 + xi);    % [W/m]   prefactor for transform

% ---- Pre-allocate outputs -----------------------------------------------
T           = T_s * ones(N, N, Nz, 'single');
u_kirchhoff = zeros(N, N, Nz, 'single');

% ---- Crystal aperture (pixels inside crystal cross-section) -------------
in_crystal = (R_g <= R_c);

% ---- FDM coefficient matrices for the 2-D Poisson equation -------------
% ∇²(Δu) = −Q  on the square (x,y) grid with dx spacing.
% We use a 5-point Laplacian stencil on the N×N grid.
% Boundary conditions on the perimeter of the grid: natural (du/dn = 0)
% at the grid boundary (>> R_c), plus Robin at the crystal surface.
%
% For efficiency, the same sparse matrix A is reused for all z-slices.
n2   = N * N;
% Indices: pixel (i,j) → linear index k = (j-1)*N + i
ii   = reshape(1:n2, N, N);

% Diagonal (centre pixel)
d_vals = -4 * ones(n2, 1);

% Off-diagonal connections (interior)
% right (+x): (i, j) → (i+1, j)
% left  (−x): (i, j) → (i−1, j)
% up    (+y): (i, j) → (i, j+1)
% down  (−y): (i, j) → (i, j−1)

% Build sparse Laplacian via diagonals
e  = ones(n2,1);
% Shift in x (periodic wrap handled by masking)
Ax = spdiags([e, e], [-1, 1], n2, n2);
% Fix wrap-around between columns: (i=N, j) and (i=1, j+1) → NOT neighbours
% These are at linear positions j*N and j*N+1 → offset ±1 → need to zero them
rows_wrap = (N : N : n2-1)';   % last pixel of each column (except last)
Ax(sub2ind([n2,n2], rows_wrap,   rows_wrap+1)) = 0;
Ax(sub2ind([n2,n2], rows_wrap+1, rows_wrap  )) = 0;

% Shift in y
Ay = spdiags([e, e], [-N, N], n2, n2);

% Laplacian matrix (scaled by 1/dx²)
L_op = (Ax + Ay + spdiags(d_vals, 0, n2, n2)) / dx^2;

% ---- Apply Robin BC at crystal surface ----------------------------------
% For pixels adjacent to the crystal surface: add correction to the
% diagonal to enforce  k·dT/dr = −h_c·(T − T_sink)
% In terms of Δu:  d(Δu)/dr|_{surface} = −h_c·(T − T_sink)
% For small temperature rises, T − T_sink ≈ Δu/k_avg, so:
%   d(Δu)/dr|_{surface} ≈ −h_c/k_avg·Δu
% We use k_avg ≈ k0 (first-order approximation):
k_avg = k0;
h_scaled = h_c / k_avg;   % effective Robin coefficient [m⁻¹]

% Find pixels just INSIDE the crystal near the surface
% Approximate: identify pixels where any neighbour is outside the crystal.
% Implemented without Image Processing Toolbox: erode by checking 4 neighbours.
in_c   = double(in_crystal);
eroded = (circshift(in_c,[ 1,0]) & circshift(in_c,[-1,0]) & ...
          circshift(in_c,[0, 1]) & circshift(in_c,[0,-1]) & in_crystal);
surface_mask = in_crystal & ~eroded;
surf_idx     = find(surface_mask);

% Add Robin penalty to diagonal: L(k,k) -= h_scaled / dx  (1-D approximation)
L_diag_corr = sparse(surf_idx, surf_idx, -h_scaled / dx * ones(numel(surf_idx),1), n2, n2);
L_total = L_op + L_diag_corr;

% Set rows outside the crystal to identity (Δu = 0 outside)
out_idx = find(~in_crystal(:));
for k_idx = out_idx'
    L_total(k_idx, :) = 0;
    L_total(k_idx, k_idx) = 1;
end

% Factorise once (LU) – reused for all z-slices
[L_lu, U_lu, P_lu, Q_lu] = lu(L_total);

% ---- Solve slice by slice -----------------------------------------------
for kz = 1 : Nz
    % Source term: RHS = Q/dx² ... no, Poisson: ∇²u = −Q → L·u = −Q
    Q_slice = double(Q_vol_3d(:,:,kz));   % [N×N]  [W/m³]

    % Build RHS vector: −Q inside crystal, 0 outside
    rhs = zeros(n2, 1);
    rhs(in_crystal(:)) = -Q_slice(in_crystal(:));

    % Solve the linear system
    du_vec = Q_lu * (U_lu \ (L_lu \ (P_lu * rhs)));
    du     = reshape(du_vec, N, N);     % Δu [W/m]

    % Inverse Kirchhoff transform → T
    % u(T) = coeff · [ (T/T0)^(1+ξ) − (T_sink/T0)^(1+ξ) ]
    % Since u=0 at T=T_sink, du is the Kirchhoff variable directly.
    % T = T0 · [ du/coeff + (T_sink/T0)^(1+ξ) ]^{1/(1+ξ)}
    base   = (T_s / T0)^(1 + xi);
    arg    = max(du ./ coeff + base, base);   % floor at T_sink (arg ≥ base)
    T_slice = T0 .* arg .^ (1/(1+xi));
    T_slice(~in_crystal) = T_s;               % outside crystal: T = T_sink

    T(:,:,kz)           = single(T_slice);
    u_kirchhoff(:,:,kz) = single(du);
end

end
